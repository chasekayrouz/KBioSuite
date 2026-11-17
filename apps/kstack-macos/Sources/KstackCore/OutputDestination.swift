import Foundation

/// A new stack follows its first input until the user chooses a destination.
public struct OutputDestination {
    public private(set) var url: URL?
    public private(set) var isExplicit = false
    public init() {}

    public mutating func follow(firstInput: URL?) {
        guard !isExplicit else { return }
        guard let firstInput else { url = nil; return }
        let folder = firstInput.deletingLastPathComponent()
        var candidate = folder.appendingPathComponent("Kstack.html")
        var suffix = 2
        while FileManager.default.fileExists(atPath: candidate.path) {
            candidate = folder.appendingPathComponent("Kstack \(suffix).html")
            suffix += 1
        }
        url = candidate
    }

    public mutating func choose(_ url: URL?) {
        self.url = url
        isExplicit = url != nil
    }
}
