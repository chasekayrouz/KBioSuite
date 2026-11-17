import SwiftUI
import WebKit

struct WebEditor: NSViewRepresentable {
    @ObservedObject var model: AppModel
    func makeCoordinator() -> Coordinator { Coordinator(model) }
    func makeNSView(context: Context) -> WKWebView {
        let configuration = WKWebViewConfiguration()
        configuration.userContentController.add(context.coordinator, name: "kstack")
        configuration.userContentController.addUserScript(WKUserScript(source: "window.addEventListener('DOMContentLoaded',()=>window.Kstack?.setNativeMode(true));", injectionTime: .atDocumentEnd, forMainFrameOnly: true))
        let view = WKWebView(frame: .zero, configuration: configuration)
        view.navigationDelegate = context.coordinator
        view.setValue(false, forKey: "drawsBackground")
        model.webView = view
        return view
    }
    func updateNSView(_ view: WKWebView, context: Context) {
        if context.coordinator.loadedHTML != model.html {
            context.coordinator.loadedHTML = model.html
            view.loadHTMLString(model.html, baseURL: nil)
        }
    }
    static func dismantleNSView(_ nsView: WKWebView, coordinator: Coordinator) {
        nsView.configuration.userContentController.removeScriptMessageHandler(forName: "kstack")
    }
    final class Coordinator: NSObject, WKScriptMessageHandler, WKNavigationDelegate {
        let model: AppModel
        var loadedHTML = ""
        init(_ model: AppModel) { self.model = model }
        func userContentController(_ userContentController: WKUserContentController, didReceive message: WKScriptMessage) {
            guard message.frameInfo.isMainFrame, let value = message.body as? [String: Any] else { return }
            model.receive(value)
        }
        func webView(_ webView: WKWebView, didFinish navigation: WKNavigation!) {
            webView.evaluateJavaScript("window.Kstack?.setNativeMode(true)")
        }
        func webView(_ webView: WKWebView, decidePolicyFor navigationAction: WKNavigationAction, decisionHandler: @escaping (WKNavigationActionPolicy) -> Void) {
            let scheme = navigationAction.request.url?.scheme
            decisionHandler(scheme == "about" || scheme == nil ? .allow : .cancel)
        }
    }
}
