// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "Kstack",
    platforms: [.macOS("26.0")],  // kbiosuite.version.MINIMUM_MACOS; checked by check_versions.py
    products: [.executable(name: "Kstack", targets: ["Kstack"])],
    targets: [
        .target(name: "KstackCore"),
        .executableTarget(name: "Kstack", dependencies: ["KstackCore"]),
        .testTarget(name: "KstackCoreTests", dependencies: ["KstackCore"])
    ],
    swiftLanguageModes: [.v5]
)
