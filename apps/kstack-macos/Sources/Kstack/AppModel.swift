import AppKit
import Combine
import UniformTypeIdentifiers
import WebKit
import KstackCore

struct ImportEntry: Identifiable {
    var path: String
    var name: String
    var cluster: JSONValue?
    var error: String?
    var circularChoice: String = ""
    var id: String { path }
    var circular: Bool { cluster?["topology"].string == "circular" }
    var omitted: Bool { circular && circularChoice == "omit" }
    var warningCount: Int {
        (cluster?["warnings"].array.count ?? 0) + (cluster?["cds"].array.reduce(0) { $0 + $1["warnings"].array.count } ?? 0)
    }
}

@MainActor
final class AppModel: ObservableObject {
    @Published var document = DocumentState()
    @Published var inputs: [ImportEntry] = [] {
        didSet { destination.follow(firstInput: inputs.first.map { URL(fileURLWithPath: $0.path) }) }
    }
    @Published var html = ""
    @Published var busy = false
    @Published var smokeRunning = false
    @Published var status = ""
    @Published var errorMessage: String?
    @Published var recent: [RecentEntry] = []
    @Published var identity = 30.0
    @Published var coverage = 70.0
    @Published var preserve = false
    @Published var connectors = true
    @Published private var destination = OutputDestination()
    var outputURL: URL? {
        get { destination.url }
        set { destination.choose(newValue) }
    }
    @Published var recomputeIdentity = 30.0
    @Published var recomputeCoverage = 70.0
    private(set) var backend: Backend?
    let recents = RecentStore()
    weak var webView: WKWebView?
    weak var window: NSWindow?
    let shutdown = ShutdownCoordinator()
    var smokeMode = false
    var onReady: (() -> Void)?

    init() {
        recent = recents.entries
        do { backend = try Backend.resolve() }
        catch { errorMessage = error.localizedDescription }
    }
    var canGenerate: Bool {
        backend != nil && !busy && !inputs.isEmpty && inputs.count <= 50 && outputURL != nil &&
        inputs.allSatisfy { $0.error == nil && (!$0.circular || !$0.circularChoice.isEmpty) } &&
        inputs.contains { !$0.omitted } && identity >= 0 && identity <= 100 && coverage >= 0 && coverage <= 100
    }
    func request(_ object: [String: JSONValue]) async throws -> JSONValue {
        guard let backend else { throw BackendError("Python/Biopython could not be located.") }
        return try await backend.request(.object(object))
    }
    func report(_ error: Error) { errorMessage = error.localizedDescription }
    func addFiles() {
        guard !busy else { return }
        let panel = NSOpenPanel(); panel.allowsMultipleSelection = true; panel.canChooseDirectories = false
        panel.title = "Add GenBank Files"; panel.prompt = "Add Files"
        panel.allowedContentTypes = ["gb", "gbk", "gbff", "genbank"].compactMap { UTType(filenameExtension: $0) }
        if panel.runModal() == .OK { Task { await importURLs(panel.urls) } }
    }
    func importURLs(_ urls: [URL]) async {
        guard !busy else { return }
        let existing = Set(inputs.map(\.path))
        var added = Set<String>()
        let paths = urls.map { $0.standardizedFileURL.path }.filter { !existing.contains($0) && added.insert($0).inserted }
        guard !paths.isEmpty else { return }
        guard inputs.count + paths.count <= 50 else { errorMessage = "Kstack supports at most 50 files per stack."; return }
        let old = Dictionary(uniqueKeysWithValues: inputs.map { ($0.path, $0) })
        let allPaths = inputs.map(\.path) + paths
        busy = true; status = "Reading GenBank features and proteins…"; defer { busy = false }
        do {
            let result = try await request(["action": .string("inspect"), "paths": .array(allPaths.map(JSONValue.string))])
            inputs = result["inputs"].array.map { value in
                let path = value["path"].string ?? ""
                let cluster = value["cluster"] == .null ? nil : value["cluster"]
                let previous = old[path]
                return ImportEntry(path: path, name: previous?.name ?? cluster?["name"].string ?? URL(fileURLWithPath: path).deletingPathExtension().lastPathComponent,
                                   cluster: cluster, error: value["error"].string, circularChoice: previous?.circularChoice ?? "")
            }
            status = "\(inputs.count) files inspected"
        } catch { report(error) }
    }
    func chooseOutput() {
        let panel = NSSavePanel(); panel.title = "Save New Stack"; panel.allowedContentTypes = [.html]
        panel.nameFieldStringValue = outputURL?.lastPathComponent ?? "Kstack.html"
        if let outputURL { panel.directoryURL = outputURL.deletingLastPathComponent() }
        if panel.runModal() == .OK { outputURL = panel.url }
    }
    func moveInput(_ id: String, by delta: Int) {
        guard let index = inputs.firstIndex(where: { $0.id == id }) else { return }
        let next = index + delta
        guard inputs.indices.contains(next) else { return }; inputs.swapAt(index, next)
    }
    func generate() async {
        guard canGenerate, let outputURL else { return }
        busy = true; status = "Finding homologs with MMseqs2…"; defer { busy = false }
        do {
            let clusters = inputs.filter { !$0.omitted }.compactMap { entry -> JSONValue? in
                guard var cluster = entry.cluster else { return nil }
                let name = entry.name.trimmingCharacters(in: .whitespacesAndNewlines)
                cluster["name"] = .string(name.isEmpty ? URL(fileURLWithPath: entry.path).deletingPathExtension().lastPathComponent : name)
                if entry.circular { cluster["linearized_at_origin"] = .bool(true) }
                return cluster
            }
            let result = try await request(["action": .string("generate"), "clusters": .array(clusters),
                "title": .string(outputURL.deletingPathExtension().lastPathComponent),
                "settings": .object(["identity": .number(identity/100), "coverage": .number(coverage/100)]),
                "preserve": .bool(preserve), "connectors": .bool(connectors)])
            guard let html = result["html"].string else { throw BackendError("Missing HTML from backend.") }
            if !destination.isExplicit && FileManager.default.fileExists(atPath: outputURL.path) {
                throw BackendError("The default output already exists. Choose a destination before generating this stack.")
            }
            try Data(html.utf8).write(to: outputURL, options: .atomic)
            document.open(project: result["project"], at: outputURL)
            self.html = html; setAnalysisControls(); remember(outputURL); updateWindow()
            status = "Generated in \(String(format: "%0.2f", result["elapsed_seconds"].number ?? 0)) seconds"
        } catch { report(error) }
    }
    func browseOpen() {
        guard !busy else { return }
        let panel = NSOpenPanel(); panel.allowedContentTypes = [.html]; panel.title = "Open Kstack HTML"
        if panel.runModal() == .OK, let url = panel.url { Task { await open(url) } }
    }
    func open(_ url: URL) async {
        guard !busy else { return }
        guard await protectUnsaved() else { return }
        busy = true; status = "Opening stack…"; defer { busy = false }
        do {
            let result = try await request(["action": .string("open"), "path": .string(url.path)])
            // Never execute JavaScript supplied by an arbitrary HTML project. Read data,
            // validate it, and regenerate a renderer from this app's trusted resources.
            let rendered = try await request(["action": .string("serialize"), "project": result["project"]])
            document.open(project: result["project"], at: url); html = rendered["html"].string ?? ""
            inputs = []; setAnalysisControls(); remember(url); updateWindow(); status = "Opened \(url.lastPathComponent)"
        } catch { report(error) }
    }
    func setAnalysisControls() {
        recomputeIdentity = (document.project?["analysis"]["identity"].number ?? 0.3)*100
        recomputeCoverage = (document.project?["analysis"]["coverage"].number ?? 0.7)*100
    }
    func remember(_ url: URL) {
        if !smokeMode { recents.remember(url); recent = recents.entries }
    }
    func removeRecent(_ entry: RecentEntry) { recents.remove(entry); recent = recents.entries }
    func syncWebProject() async throws {
        if let webView, document.project != nil {
            let value = try await webView.evaluateJavaScript("window.Kstack.getProject()")
            if let value { document.update(try JSONValue.from(value)) }
        }
    }
    @discardableResult
    func save(as newLocation: Bool = false) async -> Bool {
        guard !busy, document.project != nil else { return false }
        var url = document.fileURL
        if newLocation || url == nil {
            let panel = NSSavePanel(); panel.allowedContentTypes = [.html]; panel.title = "Save Stack As"
            panel.nameFieldStringValue = document.fileURL?.lastPathComponent ?? "Kstack.html"
            panel.directoryURL = document.fileURL?.deletingLastPathComponent()
            guard panel.runModal() == .OK else { return false }; url = panel.url
        }
        guard let destination = url else { return false }
        busy = true; status = "Saving HTML project…"; defer { busy = false }
        do {
            try await syncWebProject(); guard let snapshot = document.project else { return false }
            let result = try await request(["action": .string("serialize"), "project": snapshot])
            guard let html = result["html"].string else { throw BackendError("Missing saved HTML.") }
            try Data(html.utf8).write(to: destination, options: .atomic)
            document.saved(at: destination, project: snapshot); remember(destination); updateWindow(); status = "Saved \(destination.lastPathComponent)"
            return true
        } catch { report(error); return false }
    }
    func protectUnsaved() async -> Bool {
        guard document.isDirty else { return true }
        let alert = NSAlert(); alert.messageText = "Save changes to “\(document.title)”?"
        alert.informativeText = "Your changes will be lost if you don’t save them."
        alert.addButton(withTitle: "Save"); alert.addButton(withTitle: "Cancel"); alert.addButton(withTitle: "Don’t Save")
        let decision: SaveDecision
        switch alert.runModal() {
        case .alertFirstButtonReturn: decision = .save
        case .alertThirdButtonReturn: decision = .discard
        default: decision = .cancel
        }
        return await decision.resolve { await self.save() }
    }
    func newStack() async {
        guard !busy, await protectUnsaved() else { return }
        document.clear(); html = ""; inputs = []; outputURL = nil; status = ""; updateWindow()
    }
    func closeStack() async {
        guard !busy, await protectUnsaved() else { return }
        document.clear(); html = ""; inputs = []; updateWindow()
    }
    func recompute() async {
        guard !busy, document.project != nil, 0...100 ~= recomputeIdentity, 0...100 ~= recomputeCoverage else { return }
        busy = true; status = "Recomputing homologs from stored proteins…"; defer { busy = false }
        do {
            try await syncWebProject(); guard let snapshot = document.project else { return }
            let result = try await request(["action": .string("recompute"), "project": snapshot,
                "settings": .object(["identity": .number(recomputeIdentity/100), "coverage": .number(recomputeCoverage/100)])])
            let rendered = try await request(["action": .string("serialize"), "project": result["project"]])
            document.update(result["project"]); html = rendered["html"].string ?? ""; updateWindow(); status = "Homologs updated; save to keep these changes"
        } catch { report(error) }
    }
    func receive(_ message: [String: Any]) {
        switch message["action"] as? String {
        case "changed":
            if let value = message["project"], let json = try? JSONValue.from(value) { document.update(json); updateWindow() }
        case "copy":
            if let text = message["content"] as? String { NSPasteboard.general.clearContents(); NSPasteboard.general.setString(text, forType: .string) }
        case "export":
            if let content = message["content"] as? String, let kind = message["kind"] as? String {
                exportText(content, kind: kind, name: message["filename"] as? String ?? "Kstack.\(kind)")
            }
        case "ready": onReady?()
        default: break
        }
    }
    func exportText(_ content: String, kind: String, name: String) {
        guard ["faa", "svg"].contains(kind) else { return }
        let panel = NSSavePanel(); panel.title = kind == "svg" ? "Export Figure as SVG" : "Export Selected Proteins"
        panel.allowedContentTypes = [UTType(filenameExtension: kind) ?? .plainText]; panel.nameFieldStringValue = name
        if panel.runModal() == .OK, let url = panel.url {
            do { try Data(content.utf8).write(to: url, options: .atomic); status = "Exported \(url.lastPathComponent)" } catch { report(error) }
        }
    }
    func exportFromMenu(_ kind: String) {
        guard !busy else { return }
        Task {
            do {
                if kind == "svg" {
                    _ = try await webView?.evaluateJavaScript("Kstack.openExportPreview(); true")
                    return
                }
                let value = try await webView?.evaluateJavaScript("Kstack.exportFASTA()")
                if let content = value as? String, !content.isEmpty { exportText(content, kind: kind, name: "Kstack.\(kind)") }
                else { errorMessage = "Select at least one protein before exporting FASTA." }
            } catch { report(error) }
        }
    }
    func updateWindow() {
        window?.isDocumentEdited = document.isDirty
        window?.title = document.project == nil ? "Kstack" : document.title
        window?.representedURL = document.fileURL
    }
}
