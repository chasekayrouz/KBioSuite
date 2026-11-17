import Foundation

public struct BackendError: LocalizedError {
    public let message: String
    public var errorDescription: String? { message }
    public init(_ message: String) { self.message = message }
}

public struct Backend {
    public let resources: URL
    public let python: URL
    public init(resources: URL, python: URL) { self.resources = resources; self.python = python }

    public static func resolve(environment: [String: String] = ProcessInfo.processInfo.environment,
                               bundle: Bundle = .main) throws -> Backend {
        let fm = FileManager.default
        var roots: [URL] = []
        if let path = environment["KSTACK_RESOURCES"] { roots.append(URL(fileURLWithPath: path)) }
        if let resources = bundle.resourceURL { roots.append(resources) }
        var ancestor = URL(fileURLWithPath: CommandLine.arguments[0]).standardizedFileURL.deletingLastPathComponent()
        for _ in 0..<10 { roots.append(ancestor); ancestor.deleteLastPathComponent() }
        roots.append(URL(fileURLWithPath: fm.currentDirectoryPath))
        for root in roots {
            let frozen = root.appendingPathComponent("runtime/KBioRuntime")
            if fm.isExecutableFile(atPath: frozen.path) {
                return Backend(resources: root, python: frozen)
            }
        }
        guard let root = roots.first(where: { fm.fileExists(atPath: $0.appendingPathComponent("src/kstack/cli.py").path) }) else {
            throw BackendError("Backend resources not found. Run scripts/build-install.sh.")
        }
        var candidates: [String] = []
        if let explicit = environment["KSTACK_PYTHON"] { candidates.append(explicit) }
        candidates.append(root.appendingPathComponent(".venv/bin/python").path)
        for dir in (environment["PATH"] ?? "/usr/bin:/bin").split(separator: ":") {
            candidates.append(String(dir) + "/python3")
        }
        candidates += ["/opt/homebrew/bin/python3", "/usr/local/bin/python3", "/usr/bin/python3"]
        for candidate in candidates {
            guard fm.isExecutableFile(atPath: candidate) else { continue }
            let check = Process(); check.executableURL = URL(fileURLWithPath: candidate)
            check.arguments = ["-B", "-c", "import Bio; import sys; assert sys.version_info >= (3,12)"]
            check.standardOutput = FileHandle.nullDevice; check.standardError = FileHandle.nullDevice
            do { try check.run(); check.waitUntilExit(); if check.terminationStatus == 0 { return Backend(resources: root, python: URL(fileURLWithPath: candidate)) } } catch { continue }
        }
        throw BackendError("Python 3.12+ with Biopython is required. Run uv sync or set KSTACK_PYTHON.")
    }

    public func request(_ value: JSONValue) async throws -> JSONValue {
        try await Task.detached(priority: .userInitiated) { try self.run(value) }.value
    }
    public func run(_ value: JSONValue) throws -> JSONValue {
        let process = Process(), input = Pipe(), output = Pipe()
        let temporary = FileManager.default.temporaryDirectory.appendingPathComponent("kstack-stderr-\(UUID().uuidString)")
        FileManager.default.createFile(atPath: temporary.path, contents: nil)
        let error = try FileHandle(forWritingTo: temporary)
        defer { try? error.close(); try? FileManager.default.removeItem(at: temporary) }
        process.executableURL = python
        let frozen = python.lastPathComponent == "KBioRuntime"
        process.arguments = frozen ? [] : ["-B", "-m", "kstack.cli"]
        var environment = ProcessInfo.processInfo.environment
        // Signed app resources are immutable; all scientific scratch files remain temporary.
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        environment["PYTHONUNBUFFERED"] = "1"
        if frozen {
            environment.removeValue(forKey: "PYTHONPATH")
            environment["KSTACK_RESOURCES"] = resources.appendingPathComponent("kstack").path
            environment["KSTACK_MMSEQS"] = resources.appendingPathComponent("runtime/mmseqs/bin/mmseqs").path
        } else {
            environment["PYTHONPATH"] = resources.appendingPathComponent("src").path
            environment["KSTACK_RESOURCES"] = resources.appendingPathComponent("src/kstack/resources").path
        }
        process.environment = environment
        process.standardInput = input; process.standardOutput = output; process.standardError = error
        try process.run()
        let response = ResponseBuffer()
        let group = DispatchGroup(); group.enter()
        DispatchQueue.global(qos: .userInitiated).async {
            response.data = output.fileHandleForReading.readDataToEndOfFile(); group.leave()
        }
        do { try input.fileHandleForWriting.write(contentsOf: value.encoded()); try input.fileHandleForWriting.close() }
        catch { process.terminate(); throw error }
        process.waitUntilExit(); group.wait()
        let envelope = try JSONValue.decode(response.data)
        guard envelope["ok"].bool == true else {
            throw BackendError(envelope["error"]["message"].string ?? "Biological backend failed.")
        }
        return envelope["result"]
    }
}
private final class ResponseBuffer: @unchecked Sendable { var data = Data() }
