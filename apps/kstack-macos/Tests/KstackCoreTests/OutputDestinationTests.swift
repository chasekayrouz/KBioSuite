import Foundation
import XCTest
@testable import KstackCore

final class OutputDestinationTests: XCTestCase {
    func testFollowsFirstInputAndClearsWithoutAnOverride() {
        var destination = OutputDestination()
        XCTAssertNil(destination.url)
        destination.follow(firstInput: URL(fileURLWithPath: "/first folder/a.gbk"))
        XCTAssertEqual(destination.url?.path, "/first folder/Kstack.html")
        destination.follow(firstInput: URL(fileURLWithPath: "/second/b.gb"))
        XCTAssertEqual(destination.url?.path, "/second/Kstack.html")
        destination.follow(firstInput: nil)
        XCTAssertNil(destination.url)
    }

    func testExplicitChoiceOverridesDefault() {
        var destination = OutputDestination()
        let input = URL(fileURLWithPath: "/first folder/a.gbk")
        destination.follow(firstInput: input)
        destination.choose(destination.url)
        destination.follow(firstInput: URL(fileURLWithPath: "/second/b.gb"))
        XCTAssertEqual(destination.url?.path, "/first folder/Kstack.html")
        destination.follow(firstInput: nil)
        XCTAssertNotNil(destination.url)
        destination.choose(nil)
        destination.follow(firstInput: input)
        XCTAssertFalse(destination.isExplicit)
        XCTAssertEqual(destination.url?.path, "/first folder/Kstack.html")
    }

    func testDefaultAvoidsExistingFilesWithoutChangingThem() throws {
        let folder = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: folder, withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(at: folder) }
        for name in ["Kstack.html", "Kstack 2.html"] {
            try Data("existing stack".utf8).write(to: folder.appendingPathComponent(name))
        }
        var destination = OutputDestination()
        destination.follow(firstInput: folder.appendingPathComponent("input.gbk"))
        XCTAssertEqual(destination.url, folder.appendingPathComponent("Kstack 3.html"))
        XCTAssertEqual(try String(contentsOf: folder.appendingPathComponent("Kstack.html"), encoding: .utf8), "existing stack")
    }
}
