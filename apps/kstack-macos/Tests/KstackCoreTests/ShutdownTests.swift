import XCTest
@testable import KstackCore

final class ShutdownTests: XCTestCase {
    @MainActor
    func testCleanOrManuallySavedSessionNeedsNoDecision() {
        var document = DocumentState()
        let shutdown = ShutdownCoordinator()
        XCTAssertEqual(shutdown.request(isDirty: document.isDirty), .allow)
        document.create(project: .object([:]))
        document.saved(at: URL(fileURLWithPath: "/tmp/test.html"), project: .object([:]))
        XCTAssertEqual(shutdown.request(isDirty: document.isDirty), .allow)
    }

    @MainActor
    func testWindowCloseAndDirectQuitShareOnePendingDecision() {
        let shutdown = ShutdownCoordinator()
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
        // The same transaction covers close, quit and a reentrant quit during Save.
        XCTAssertEqual(shutdown.request(isDirty: true), .pending)
        XCTAssertEqual(shutdown.request(isDirty: true, blocked: true), .pending)
        shutdown.complete(allowed: true)
        // Don't Save deliberately leaves the document dirty, but has authorized exit.
        XCTAssertEqual(shutdown.request(isDirty: true), .allow)
    }

    @MainActor
    func testSaveRunsExactlyOnceThenAllowsExit() async {
        let shutdown = ShutdownCoordinator()
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
        var saves = 0
        let allowed = await SaveDecision.save.resolve {
            saves += 1
            XCTAssertEqual(shutdown.request(isDirty: true, blocked: true), .pending)
            return true
        }
        shutdown.complete(allowed: allowed)
        XCTAssertEqual(saves, 1)
        XCTAssertEqual(shutdown.request(isDirty: false), .allow)
    }

    @MainActor
    func testDiscardNeverSavesOrAsksAgain() async {
        let shutdown = ShutdownCoordinator()
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
        let allowed = await SaveDecision.discard.resolve { XCTFail("Unexpected save"); return false }
        shutdown.complete(allowed: allowed)
        XCTAssertTrue(allowed)
        XCTAssertEqual(shutdown.request(isDirty: true), .allow)
    }

    @MainActor
    func testCancelKeepsProtectionForTheNextShutdown() async {
        let shutdown = ShutdownCoordinator()
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
        let allowed = await SaveDecision.cancel.resolve { XCTFail("Unexpected save"); return true }
        shutdown.complete(allowed: allowed)
        XCTAssertFalse(allowed)
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
    }

    @MainActor
    func testFailedOrCancelledSaveAbortsShutdown() async {
        let shutdown = ShutdownCoordinator()
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
        let allowed = await SaveDecision.save.resolve { false }
        shutdown.complete(allowed: allowed)
        XCTAssertFalse(allowed)
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
    }

    @MainActor
    func testBusyWorkBlocksWithoutConsumingTheDecision() {
        let shutdown = ShutdownCoordinator()
        XCTAssertEqual(shutdown.request(isDirty: true, blocked: true), .cancel)
        XCTAssertEqual(shutdown.request(isDirty: true), .confirm)
    }
}
