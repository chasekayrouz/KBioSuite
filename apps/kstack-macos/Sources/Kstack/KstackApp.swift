import SwiftUI
import AppKit

@main
struct KstackApp: App {
    @NSApplicationDelegateAdaptor(AppDelegate.self) var delegate
    @StateObject private var model = AppModel()
    var body: some Scene {
        Window("Kstack", id: "main") {
            ContentView(model: model).onAppear { delegate.model = model; model.window?.makeKeyAndOrderFront(nil); runStartup() }
        }.defaultSize(width: 1320, height: 850).windowResizability(.contentMinSize)
            .commands {
                CommandGroup(replacing: .newItem) {
                    Button("New Stack") { Task { await model.newStack() } }.keyboardShortcut("n").disabled(model.busy || model.smokeRunning)
                    Button("Open Stack…") { model.browseOpen() }.keyboardShortcut("o").disabled(model.busy || model.smokeRunning)
                    Menu("Open Recent") {
                        ForEach(model.recent) { entry in
                            Button(entry.name) { Task { await model.open(entry.resolvedURL()) } }.disabled(!entry.exists || model.busy)
                        }
                    }
                }
                CommandGroup(replacing: .saveItem) {
                    Button("Save") { Task { await model.save() } }.keyboardShortcut("s").disabled(model.document.project == nil || model.busy || model.smokeRunning)
                    Button("Save As…") { Task { await model.save(as: true) } }.keyboardShortcut("s", modifiers: [.command, .shift]).disabled(model.document.project == nil || model.busy || model.smokeRunning)
                }
                CommandGroup(after: .saveItem) {
                    Divider()
                    Button("Export SVG…") { model.exportFromMenu("svg") }.disabled(model.document.project == nil || model.busy || model.smokeRunning)
                    Button("Export Protein FASTA…") { model.exportFromMenu("faa") }.disabled(model.document.project == nil || model.busy || model.smokeRunning)
                    Divider()
                    Button("Close Stack") { Task { await model.closeStack() } }.disabled(model.document.project == nil || model.busy || model.smokeRunning)
                }
            }
    }
    private func runStartup() {
        guard !delegate.started else { return }; delegate.started = true
        let args = CommandLine.arguments
        if let i = args.firstIndex(of: "--open"), args.indices.contains(i+1) {
            Task { await model.open(URL(fileURLWithPath: args[i+1])) }
        }
        if let i = args.firstIndex(of: "--smoke"), args.indices.contains(i+1) {
            Task { await NativeSmoke.run(model: model, root: URL(fileURLWithPath: args[i+1])) }
        }
    }
}

@MainActor
final class AppDelegate: NSObject, NSApplicationDelegate {
    weak var model: AppModel?
    var started = false
    func applicationDidFinishLaunching(_ notification: Notification) { NSApp.setActivationPolicy(.regular); NSApp.activate(ignoringOtherApps: true) }
    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool { true }
    func applicationShouldTerminate(_ sender: NSApplication) -> NSApplication.TerminateReply {
        guard let model else { return .terminateNow }
        switch model.shutdown.request(isDirty: model.document.isDirty, blocked: model.busy || model.smokeRunning) {
        case .allow: return .terminateNow
        case .cancel: return .terminateCancel
        case .pending: return .terminateLater
        case .confirm:
            Task {
                let confirmed = await model.protectUnsaved()
                model.shutdown.complete(allowed: confirmed)
                sender.reply(toApplicationShouldTerminate: confirmed)
            }
            return .terminateLater
        }
    }
    func application(_ application: NSApplication, open urls: [URL]) {
        if let url = urls.first { Task { await model?.open(url) } }
    }
}

struct WindowHook: NSViewRepresentable {
    @ObservedObject var model: AppModel
    func makeCoordinator() -> Coordinator { Coordinator(model) }
    func makeNSView(context: Context) -> NSView { let view = NSView(); attach(view, context.coordinator); return view }
    func updateNSView(_ view: NSView, context: Context) { attach(view, context.coordinator) }
    func attach(_ view: NSView, _ coordinator: Coordinator) {
        DispatchQueue.main.async { if let window = view.window { model.window = window; window.delegate = coordinator; model.updateWindow() } }
    }
    final class Coordinator: NSObject, NSWindowDelegate {
        let model: AppModel
        init(_ model: AppModel) { self.model = model }
        func windowShouldClose(_ sender: NSWindow) -> Bool {
            // This is a single-window application and closing it quits the app.
            // Termination owns the decision; never confirm and then close into
            // a second applicationShouldTerminate confirmation.
            NSApp.terminate(sender)
            return false
        }
    }
}
