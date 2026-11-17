import XCTest
@testable import KstackCore

final class CoreTests: XCTestCase {
    let project: JSONValue = .object(["format": .string("Kstack"),"schema_version": .number(1)])
    func testJSONRoundtrip() throws {
        let value: JSONValue = .object(["string": .string("</script> α"),"number":.number(1.25),"null":.null,"bool":.bool(false),"array":.array([.number(2)])])
        XCTAssertEqual(try JSONValue.decode(value.encoded()),value)
    }
    func testJSONMutation() {
        var p = project; p["nested"] = .object(["offset":.number(-1250)])
        XCTAssertEqual(p["nested"]["offset"].number,-1250)
    }
    func testNewDirty() {
        var state = DocumentState(); state.create(project: project)
        XCTAssertTrue(state.isDirty); XCTAssertNil(state.fileURL)
    }
    func testOpenClean() {
        var state = DocumentState(); let url = URL(fileURLWithPath:"/tmp/Stack.html")
        state.open(project:project,at:url); XCTAssertFalse(state.isDirty); XCTAssertEqual(state.title,"Stack")
    }
    func testPersistentChangesDirty() {
        var state = DocumentState(); state.open(project: project,at:URL(fileURLWithPath:"/tmp/a.html"))
        var edited=project;edited["layout"] = .string("changed");state.update(edited)
        XCTAssertTrue(state.isDirty)
    }
    func testIdenticalUpdateStaysClean() {
        var state=DocumentState();state.open(project:project,at:URL(fileURLWithPath:"/tmp/a.html"));state.update(project)
        XCTAssertFalse(state.isDirty)
    }
    func testExplicitSaveAndSaveAs() {
        var state=DocumentState();state.create(project:project)
        let url=URL(fileURLWithPath:"/tmp/NewName.html");state.saved(at:url,project:project)
        XCTAssertFalse(state.isDirty);XCTAssertEqual(state.fileURL,url)
    }
    func testClearState() {
        var state=DocumentState();state.create(project:project);state.clear()
        XCTAssertFalse(state.isDirty);XCTAssertNil(state.project)
    }
    func testRecentReferencesAndDeduplication() throws {
        let suite="Kstack.tests.\(UUID().uuidString)",defaults=UserDefaults(suiteName:suite)!
        defer { defaults.removePersistentDomain(forName:suite) }
        let store=RecentStore(defaults:defaults)
        let a=URL(fileURLWithPath:"/tmp/a.html"),b=URL(fileURLWithPath:"/tmp/b.html")
        store.remember(a);store.remember(b);store.remember(a)
        XCTAssertEqual(store.entries.map(\.path),[a.path,b.path])
        store.remove(store.entries[0]);XCTAssertEqual(store.entries.count,1)
        let json=String(decoding: defaults.data(forKey:"RecentStacks.v1")!,as:UTF8.self)
        XCTAssertFalse(json.contains("schema_version"))
    }
    func testRecentMissingFileSafe() {
        let entry=RecentEntry(path:"/tmp/kstack-nonexistent-\(UUID().uuidString).html")
        XCTAssertFalse(entry.exists);XCTAssertEqual(entry.resolvedURL().path,entry.path)
    }
    func testRecentBookmarkTracksMove() throws {
        let base=FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at:base,withIntermediateDirectories:true)
        defer { try? FileManager.default.removeItem(at:base) }
        let old=base.appendingPathComponent("Old.html"),new=base.appendingPathComponent("Moved.html")
        try Data("fixture".utf8).write(to:old)
        let entry=RecentEntry(path:old.path,bookmark:try old.bookmarkData(options:[.minimalBookmark],includingResourceValuesForKeys:nil,relativeTo:nil))
        try FileManager.default.moveItem(at:old,to:new)
        XCTAssertEqual(entry.resolvedURL().resolvingSymlinksInPath(),new.resolvingSymlinksInPath());XCTAssertTrue(entry.exists)
    }
    func testBackendIndependentJSONInterface() throws {
        let root=URL(fileURLWithPath:#filePath).deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent()
        let backend=try Backend.resolve(environment:["KSTACK_RESOURCES":root.path])
        let result=try backend.run(.object(["action":.string("inspect"),"paths":.array([.string(backend.resources.appendingPathComponent("tests/kstack/fixtures/realistic/A. alliaceus CBS 536.65.gbk").path)])]))
        XCTAssertEqual(result["inputs"].array.first?["cluster"]["cds"].array.count,9)
    }
    func testBackendLeavesResourceDirectoryImmutable() throws {
        let root=URL(fileURLWithPath:#filePath).deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent()
        let existing=try Backend.resolve(environment:["KSTACK_RESOURCES":root.path])
        let temporary=FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at:temporary,withIntermediateDirectories:true)
        defer { try? FileManager.default.removeItem(at:temporary) }
        try FileManager.default.copyItem(at:root.appendingPathComponent("src"),to:temporary.appendingPathComponent("src"))
        let cache=temporary.appendingPathComponent("src/kstack/__pycache__")
        if FileManager.default.fileExists(atPath:cache.path) { try FileManager.default.removeItem(at:cache) }
        let backend=Backend(resources:temporary,python:existing.python)
        _=try backend.run(.object(["action":.string("inspect"),"paths":.array([.string(root.appendingPathComponent("tests/kstack/fixtures/synthetic/edge-cases.gbk").path)])]))
        XCTAssertFalse(FileManager.default.fileExists(atPath:cache.path))
    }
}
