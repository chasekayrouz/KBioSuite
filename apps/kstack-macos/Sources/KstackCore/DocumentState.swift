import Foundation

public struct DocumentState {
    public private(set) var project: JSONValue?
    public private(set) var fileURL: URL?
    public private(set) var isDirty = false
    public init() {}
    public mutating func open(project: JSONValue, at url: URL) {
        self.project = project; fileURL = url; isDirty = false
    }
    public mutating func create(project: JSONValue) {
        self.project = project; fileURL = nil; isDirty = true
    }
    public mutating func update(_ project: JSONValue) {
        if self.project != project { self.project = project; isDirty = true }
    }
    public mutating func saved(at url: URL, project: JSONValue) {
        self.project = project; fileURL = url; isDirty = false
    }
    public mutating func clear() { self = DocumentState() }
    public var title: String { fileURL?.deletingPathExtension().lastPathComponent ?? "Untitled Stack" }
}

public struct RecentEntry: Codable, Identifiable, Equatable {
    public var path: String
    public var bookmark: Data?
    public var id: String { path }
    public var name: String { URL(fileURLWithPath: path).deletingPathExtension().lastPathComponent }
    public init(path: String, bookmark: Data? = nil) { self.path = path; self.bookmark = bookmark }
    public func resolvedURL() -> URL {
        if let bookmark {
            var stale = false
            if let resolved = try? URL(resolvingBookmarkData: bookmark, options: [.withoutUI], relativeTo: nil, bookmarkDataIsStale: &stale) { return resolved }
        }
        return URL(fileURLWithPath: path)
    }
    public var exists: Bool { FileManager.default.fileExists(atPath: resolvedURL().path) }
}

public final class RecentStore {
    private let defaults: UserDefaults
    private let key: String
    public var entries: [RecentEntry] {
        guard let data = defaults.data(forKey: key) else { return [] }
        return (try? JSONDecoder().decode([RecentEntry].self, from: data)) ?? []
    }
    public init(defaults: UserDefaults = .standard, key: String = "RecentStacks.v1") { self.defaults = defaults; self.key = key }
    public func remember(_ url: URL) {
        let path = url.standardizedFileURL.path
        let bookmark = try? url.bookmarkData(options: [.minimalBookmark], includingResourceValuesForKeys: nil, relativeTo: nil)
        let entry = RecentEntry(path: path, bookmark: bookmark)
        let values = [entry] + entries.filter { $0.path != path && $0.resolvedURL().standardizedFileURL.path != path }
        defaults.set(try? JSONEncoder().encode(Array(values.prefix(20))), forKey: key)
    }
    public func remove(_ entry: RecentEntry) {
        defaults.set(try? JSONEncoder().encode(entries.filter { $0.id != entry.id }), forKey: key)
    }
}
