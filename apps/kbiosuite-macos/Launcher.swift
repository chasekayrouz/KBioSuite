import SwiftUI
import AppKit

@main
struct KBioSuiteLauncher: App {
    @NSApplicationDelegateAdaptor(LauncherDelegate.self) var delegate
    var body: some Scene {
        Window("KBioSuite", id: "launcher") {
            LauncherView().onAppear { delegate.runSmokeIfRequested() }
        }.defaultSize(width: 530, height: 600).windowResizability(.contentSize)
    }
}

struct LauncherView: View {
    @State private var error: String?
    let tools = [("Kblast", "Protein to genomic contexts", "magnifyingglass"),
                 ("Kannotate", "Local Pfam and CDD annotation", "text.badge.checkmark"),
                 ("Kstack", "Gene cluster comparison and editing", "rectangle.stack"),
                 ("Kmatch", "Protein structure prediction and structural homologs", "cube.transparent")]
    var body: some View {
        VStack(alignment: .leading, spacing: 22) {
            HStack { Text("KBioSuite").font(.system(size: 29, weight: .bold)); Spacer()
                Text(Bundle.main.infoDictionary?["CFBundleShortVersionString"] as? String ?? "").foregroundStyle(.secondary) }
            Text("Choose a tool").font(.title3).foregroundStyle(.secondary)
            VStack(spacing: 12) {
                ForEach(tools, id: \.0) { name, description, symbol in
                    Button {
                        LauncherDelegate.openTool(name) { issue in error = issue }
                    } label: {
                        HStack(spacing: 17) {
                            Image(systemName: symbol).font(.system(size: 24)).frame(width: 32).foregroundStyle(Color(red: 0.15, green: 0.43, blue: 0.3))
                            VStack(alignment: .leading, spacing: 5) {
                                Text(name).font(.headline)
                                Text(description).font(.subheadline).foregroundStyle(.secondary)
                            }
                            Spacer(); Image(systemName: "arrow.up.right").foregroundStyle(.secondary)
                        }.padding(16).frame(maxWidth: .infinity)
                    }.buttonStyle(.plain).background(.background, in: RoundedRectangle(cornerRadius: 12))
                        .accessibilityLabel("Open \(name)")
                }
            }
            Button("Shared resources…") {
                let base = ProcessInfo.processInfo.environment["KBIOSUITE_RESOURCES"].map {
                    URL(fileURLWithPath: NSString(string: $0).expandingTildeInPath)
                } ?? FileManager.default.homeDirectoryForCurrentUser.appendingPathComponent("Library/Application Support/KBioSuite/Resources")
                let page = base.appendingPathComponent("index.html")
                if FileManager.default.fileExists(atPath: page.path) { NSWorkspace.shared.open(page) }
                else { error = "Shared resources have not been configured. See the KBioSuite resource guide." }
            }.accessibilityLabel("View shared bioinformatics resources")
            if let error { Text(error).foregroundStyle(.red).font(.caption) }
            Text("Each tool also opens directly from Finder or Spotlight.").font(.caption).foregroundStyle(.secondary)
        }.padding(30).frame(width: 530)
    }
}

@MainActor
final class LauncherDelegate: NSObject, NSApplicationDelegate {
    private var started = false
    /// Set when the launcher runs translocated or from a read-only volume (see InstallLocation).
    static let locationProblem = InstallLocation.problem(bundleURL: Bundle.main.bundleURL)
    func applicationDidFinishLaunching(_ notification: Notification) {
        NSApp.activate()
        guard let problem = Self.locationProblem, !CommandLine.arguments.contains("--smoke") else { return }
        DispatchQueue.main.async {
            let alert = NSAlert()
            alert.alertStyle = .warning
            alert.messageText = InstallLocation.title
            alert.informativeText = problem
            alert.addButton(withTitle: "Quit")
            alert.addButton(withTitle: "Continue")
            if alert.runModal() == .alertFirstButtonReturn { NSApp.terminate(nil) }
        }
    }
    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool { true }
    static func toolURL(_ name: String) -> URL { Bundle.main.bundleURL.deletingLastPathComponent().appendingPathComponent(name + ".app") }
    static func openTool(_ name: String, completion: @escaping @MainActor (String?) -> Void) {
        let url = toolURL(name)
        guard FileManager.default.fileExists(atPath: url.path) else {
            completion(locationProblem ?? "\(name).app is missing. Keep the five suite apps together."); return
        }
        NSWorkspace.shared.openApplication(at: url, configuration: NSWorkspace.OpenConfiguration()) { app, error in
            Task { @MainActor in completion(error?.localizedDescription ?? (app == nil ? "Could not open \(name)." : nil)) }
        }
    }
    func runSmokeIfRequested() {
        guard !started else { return }; started = true
        let args = CommandLine.arguments
        guard let i = args.firstIndex(of: "--smoke"), args.indices.contains(i+1) else { return }
        let destination = URL(fileURLWithPath: args[i+1])
        Task {
            var results: [[String: Any]] = []
            for name in ["Kblast", "Kannotate", "Kstack", "Kmatch"] {
                let error: String? = await withCheckedContinuation { continuation in
                    Self.openTool(name) { continuation.resume(returning: $0) }
                }
                results.append(["tool": name, "path": Self.toolURL(name).path, "opened": error == nil, "error": error ?? ""])
                try? await Task.sleep(for: .seconds(2))
            }
            let report: [String: Any] = ["ok": results.allSatisfy { $0["opened"] as? Bool == true }, "children": results]
            try? FileManager.default.createDirectory(at: destination.deletingLastPathComponent(), withIntermediateDirectories: true)
            try? JSONSerialization.data(withJSONObject: report, options: [.prettyPrinted, .sortedKeys]).write(to: destination)
        }
    }
}
