import Foundation
import AppKit
import WebKit
import KstackCore

/// Opt-in acceptance harness. Writes only under the supplied repository's .artifacts.
@MainActor
enum NativeSmoke {
    static func run(model: AppModel, root: URL) async {
        model.smokeMode = true
        model.smokeRunning = true
        defer { model.smokeRunning = false }
        let base = root.appendingPathComponent(".artifacts/native-smoke")
        let output = base.appendingPathComponent("run-\(UUID().uuidString)")
        var checks: [String] = []
        func check(_ condition: Bool, _ name: String) throws {
            guard condition else { throw BackendError("Native smoke failed: \(name)") }; checks.append(name)
        }
        do {
            try FileManager.default.createDirectory(at: output, withIntermediateDirectories: true)
            let files = try FileManager.default.contentsOfDirectory(at: root.appendingPathComponent("tests/kstack/fixtures/realistic"), includingPropertiesForKeys: nil).filter { $0.pathExtension == "gbk" }.sorted { $0.path < $1.path }
            await model.importURLs(Array(files.prefix(3)))
            try check(model.inputs.count == 3 && model.inputs.allSatisfy { $0.error == nil }, "Native GenBank import")
            try check(model.outputURL?.deletingLastPathComponent() == files[0].deletingLastPathComponent(), "Output defaults beside the first imported GenBank")
            try check(model.canGenerate, "Default output enables generation without choosing a folder")
            model.outputURL = output.appendingPathComponent("Stack.html")
            await model.importURLs([root.appendingPathComponent("tests/kstack/fixtures/synthetic/circular-origin.gbk"), root.appendingPathComponent("tests/kstack/fixtures/synthetic/circular-origin-minus.gbk")])
            try check(!model.canGenerate, "Circular imports require a decision")
            for i in model.inputs.indices where model.inputs[i].circular {
                model.inputs[i].circularChoice = model.inputs[i].path.contains("minus") ? "omit" : "linearize"
            }
            await model.importURLs([root.appendingPathComponent("tests/kstack/fixtures/synthetic/edge-cases.gbk")])
            try check(model.outputURL == output.appendingPathComponent("Stack.html"), "Explicit destination survives more imports")
            await model.generate()
            try check(model.document.project != nil && !model.document.isDirty, "Generate and explicitly write HTML")
            try check(model.document.project?["clusters"].array.count == 5, "Circular Linearize/Omit respected")
            try await waitForWeb(model)
            let rendered = try await model.webView?.evaluateJavaScript("document.querySelectorAll('.gene').length") as? Int ?? 0
            try check(rendered > 0, "Complete native gene figure renders")
            let malformedVisible = try await model.webView?.evaluateJavaScript("(() => { const p=Kstack.getProject(),c=p.clusters.find(c=>c.source_name==='edge-cases.gbk'); Kstack.inspectNotes(c.id); const g=c.cds.find(g=>!g.segments.length); document.querySelector(`[aria-label=\"Inspect retained CDS ${g.id}\"]`).click(); return document.getElementById('inspector-content').textContent.includes('Uninterpretable CDS location') && document.querySelector('.sequence').textContent==='No interpretable protein'; })()") as? Bool ?? false
            try check(malformedVisible, "Malformed CDS remains reachable through import notes/inspector")
            let original = try Data(contentsOf: model.outputURL!)
            guard let view = model.webView else { throw BackendError("No WKWebView") }
            let copied = try await view.evaluateJavaScript("(() => { const p=Kstack.getProject(),g=p.clusters[0].cds.find(g=>g.protein); Kstack.inspect(g.id); document.querySelector('.protein-actions button').click(); return g.protein; })()") as? String
            try check(NSPasteboard.general.string(forType: .string) == copied && !(copied ?? "").isEmpty, "One-click native protein copy")
            _ = try await view.evaluateJavaScript("(() => { const input=document.getElementById('spacing-scale'); input.value=150; input.dispatchEvent(new Event('input')); return true; })()")
            try check(model.document.isDirty && model.document.project?["visualization"]["row_spacing_scale"].number == 1.5, "Row spacing alone marks the document dirty")
            _ = try await view.evaluateJavaScript("(() => { const p=Kstack.getProject(); Kstack.selectCluster(p.clusters[0].id); Kstack.edit('reverse',p.clusters[0].id); Kstack.edit('offset',p.clusters[0].id,1250); return true; })()")
            try check(model.document.isDirty && model.window?.isDocumentEdited == true, "Structural edits mark document/window dirty")
            try check(try Data(contentsOf: model.outputURL!) == original, "No autosave")
            let aligned = try await evaluateAsync(view, "await new Promise(r=>requestAnimationFrame(r)); const p=Kstack.getProject(),reference=p.clusters.flatMap(c=>c.cds).find(g=>Object.keys(CSModel.alignmentOffsets(p,g.id)).length),offsets=CSModel.alignmentOffsets(p,reference.id); document.querySelector(`[data-gene=\"${reference.id}\"]`).dispatchEvent(new MouseEvent('contextmenu',{bubbles:true,cancelable:true,clientX:100,clientY:100})); document.querySelector('body > .gene-menu button').click(); await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))); return Object.entries(offsets).every(([id,offset])=>Kstack.getProject().layout.rows[id].offset===offset)") as? Bool
            try check(aligned == true && model.document.isDirty, "Gene menu Auto-align updates offsets before explicit Save")
            let before = model.document.project
            try check(await model.save(), "Explicit Save succeeds")
            try check(!model.document.isDirty, "Save clears dirty state")
            await model.open(model.outputURL!)
            try await waitForWeb(model)
            try check(model.document.project == before, "Reopen restores complete edited state")
            try check(!model.document.isDirty && model.document.project?["visualization"]["row_spacing_scale"].number == 1.5, "Saved row spacing reopens cleanly")
            _ = try await view.evaluateJavaScript("Kstack.fit(); Kstack.zoom(1.2); true")
            try check(!model.document.isDirty, "Camera Fit/zoom leaves the document clean")
            model.exportFromMenu("svg")
            var previewReady = false
            // Ordering, not speed, is under test; allow a busy machine up to 10 s.
            for _ in 0..<100 {
                if try await view.evaluateJavaScript("document.getElementById('export-preview').open") as? Bool == true {
                    previewReady = true; break
                }
                try await Task.sleep(nanoseconds: 100_000_000)
            }
            try check(previewReady, "Native Export SVG menu opens preview before the save panel")
            _ = try await evaluateAsync(view, "document.querySelectorAll('#export-clusters input')[1].click(); await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))); return true")
            try check(!model.document.isDirty, "Export-only visibility leaves native document clean")
            try await snapshot(view, to: output.appendingPathComponent("Export-preview-hidden-row.png"))
            _ = try await view.evaluateJavaScript("document.getElementById('export-cancel').click(); true")
            let fasta = try await model.webView?.evaluateJavaScript("Kstack.exportFASTA()") as? String ?? ""
            let svg = try await model.webView?.evaluateJavaScript("Kstack.exportSVG()") as? String ?? ""
            try check(fasta.hasPrefix(">") && fasta.contains("|"), "Native selected-cluster FASTA")
            try check(svg.contains("data-kind=\"gene\"") && !svg.contains("<g>") && !svg.contains("<image"), "Native deterministic SVG export")
            try Data(fasta.utf8).write(to: output.appendingPathComponent("Selected.faa"))
            try Data(svg.utf8).write(to: output.appendingPathComponent("Figure.svg"))
            let standalone = WKWebView(frame: NSRect(x: 0, y: 0, width: 1200, height: 800))
            let browserWindow = NSWindow(contentRect: standalone.frame, styleMask: [.titled], backing: .buffered, defer: false)
            browserWindow.contentView = standalone; browserWindow.title = "Offline HTML Acceptance"
            browserWindow.orderFront(nil)
            standalone.loadFileURL(model.outputURL!, allowingReadAccessTo: output)
            var offlineReady = false
            for _ in 0..<100 {
                if !standalone.isLoading,
                   let ready = try? await standalone.evaluateJavaScript("!!window.Kstack && document.querySelectorAll('.gene').length > 0 && !document.body.classList.contains('native')"),
                   ready as? Bool == true { offlineReady = true; break }
                try await Task.sleep(nanoseconds: 100_000_000)
            }
            try check(offlineReady, "Direct file HTML renders offline without a native bridge/server")
            let offlineSVG = try await standalone.evaluateJavaScript("Kstack.exportSVG()") as? String
            try check(offlineSVG == svg, "Standalone SVG matches native geometry")
            browserWindow.orderOut(nil); model.window?.makeKeyAndOrderFront(nil)
            let harness = try String(contentsOf: root.appendingPathComponent("tests/kstack/renderer/visual-harness.js"), encoding: .utf8)
            let geometryReport = try await evaluateAsync(view, "try { return await \(harness)() } catch(error) { return {ok:false,error:String(error)} }") as? [String: Any]
            try check(geometryReport?["ok"] as? Bool == true, "WKWebView scale/zoom/export/interaction regression harness\(geometryReport?["error"].map { ": \($0)" } ?? "")")
            for artifact in geometryReport?["exports"] as? [[String: String]] ?? [] {
                if let name = artifact["name"], let content = artifact["svg"] {
                    try Data(content.utf8).write(to: output.appendingPathComponent("Figure-\(name).svg"))
                }
            }
            var geometrySummary = geometryReport ?? [:]
            geometrySummary.removeValue(forKey: "exports")
            try JSONSerialization.data(withJSONObject: geometrySummary, options: [.prettyPrinted,.sortedKeys]).write(to: output.appendingPathComponent("geometry-report.json"))
            // Capture the actual WKWebView, including native text and HTML row
            // controls, for visual review at each representative presentation.
            _ = try await view.evaluateJavaScript("window.__smokeProject=Kstack.getProject(); true")
            for (width, height, spacing) in [(100,100,100),(50,100,100),(200,100,100),(100,50,100),(100,200,100),(65,175,100),(100,100,50),(100,100,150),(65,175,150)] {
                _ = try await evaluateAsync(view, "for(const [id,value]of [['width-scale',\(width)],['height-scale',\(height)],['spacing-scale',\(spacing)]]){const input=document.getElementById(id);input.value=value;input.dispatchEvent(new Event('input'));} const labels=document.getElementById('labels');labels.checked=true;labels.dispatchEvent(new Event('change'));await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));Kstack.fit();await new Promise(r=>requestAnimationFrame(r));return true")
                try await snapshot(view, to: output.appendingPathComponent("Live-\(width)x\(height)-spacing\(spacing).png"))
            }
            // Exercise actual WebKit page magnification as well as the figure camera.
            let presentationSVG = try await view.evaluateJavaScript("Kstack.exportSVG()") as? String
            view.pageZoom = 2
            _ = try await evaluateAsync(view, "Kstack.fit(); await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))); return true")
            let zoomSVG = try await view.evaluateJavaScript("Kstack.exportSVG()") as? String
            try check(zoomSVG == presentationSVG, "200% WebKit page zoom preserves anisotropic vector export")
            try await snapshot(view, to: output.appendingPathComponent("page-zoom-200.png"))
            view.pageZoom = 1
            _ = try await evaluateAsync(view, "Kstack.setProject(window.__smokeProject);delete window.__smokeProject;await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));return true")
            try check(try await view.evaluateJavaScript("Kstack.exportSVG()") as? String == svg, "Visual review restores the saved presentation")
            model.recomputeIdentity = 40
            await model.recompute()
            try check(model.document.project?["analysis"]["identity"].number == 0.4 && model.document.isDirty, "Native recompute from stored proteins")
            try await waitForWeb(model)
            try check(await model.save(), "Save recomputed project")
            let metrics = try await model.webView?.evaluateJavaScript("Kstack.metrics")
            let report: [String: Any] = ["ok": true, "checks": checks, "metrics": metrics ?? [:], "artifact_directory": output.path]
            try JSONSerialization.data(withJSONObject: report, options: [.prettyPrinted,.sortedKeys]).write(to: output.appendingPathComponent("report.json"))
            try JSONSerialization.data(withJSONObject: report, options: [.prettyPrinted,.sortedKeys]).write(to: base.appendingPathComponent("report.json"))
            model.status = "Native acceptance passed: \(checks.count) checks"
        } catch {
            let report: [String: Any] = ["ok": false,"checks": checks,"error":error.localizedDescription,"artifact_directory":output.path]
            try? JSONSerialization.data(withJSONObject: report, options: [.prettyPrinted,.sortedKeys]).write(to: output.appendingPathComponent("report.json"))
            try? JSONSerialization.data(withJSONObject: report, options: [.prettyPrinted,.sortedKeys]).write(to: base.appendingPathComponent("report.json"))
            model.errorMessage = error.localizedDescription
        }
    }
    private static func snapshot(_ view: WKWebView, to url: URL) async throws {
        let image = try await view.takeSnapshot(configuration: nil)
        guard let tiff = image.tiffRepresentation, let bitmap = NSBitmapImageRep(data: tiff), let png = bitmap.representation(using: .png, properties: [:]) else {
            throw BackendError("WKWebView snapshot produced no PNG")
        }
        try png.write(to: url)
    }
    private static func evaluateAsync(_ view: WKWebView, _ body: String) async throws -> Any {
        try await withCheckedThrowingContinuation { continuation in
            view.callAsyncJavaScript(body, arguments: [:], in: nil, in: .page) { result in
                continuation.resume(with: result)
            }
        }
    }
    private static func waitForWeb(_ model: AppModel) async throws {
        model.window?.makeKeyAndOrderFront(nil)
        NSApp.activate(ignoringOtherApps: true)
        for _ in 0..<100 {
            let identity = model.document.project?["analysis"]["identity"].number ?? -1
            if let view = model.webView, !view.isLoading,
               let ready = try? await view.evaluateJavaScript("!!window.Kstack && Kstack.getProject().analysis.identity === \(identity) && Kstack.metrics.frames > 0"),
               ready as? Bool == true {
                return
            }
            try await Task.sleep(nanoseconds: 100_000_000)
        }
        let diagnostic = try? await model.webView?.evaluateJavaScript("JSON.stringify({visibility:document.visibilityState,ready:document.readyState,identity:window.Kstack?.getProject().analysis.identity,metrics:window.Kstack?.metrics,genes:document.querySelectorAll('.gene').length,viewport:[innerWidth,innerHeight]})")
        throw BackendError("WKWebView did not load the project within 10 seconds: \(diagnostic as? String ?? "unavailable")")
    }
}
