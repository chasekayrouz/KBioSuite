/// A shutdown transaction shared by window-close and direct application quit.
/// Reentrant requests join the pending decision, including during an async save.
@MainActor
public final class ShutdownCoordinator {
    public enum Request: Equatable { case allow, cancel, confirm, pending }
    private enum Phase { case idle, deciding, approved }
    private var phase = Phase.idle
    public init() {}

    public func request(isDirty: Bool, blocked: Bool = false) -> Request {
        switch phase {
        case .approved: return .allow
        case .deciding: return .pending
        case .idle:
            if blocked { return .cancel }
            if !isDirty { return .allow }
            phase = .deciding
            return .confirm
        }
    }

    public func complete(allowed: Bool) {
        precondition(phase == .deciding)
        phase = allowed ? .approved : .idle
    }
}

public enum SaveDecision {
    case save, discard, cancel

    @MainActor
    public func resolve(save: () async -> Bool) async -> Bool {
        switch self {
        case .save: return await save()
        case .discard: return true
        case .cancel: return false
        }
    }
}
