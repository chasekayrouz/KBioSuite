import Foundation

/// Whether the launcher runs from a place where its sibling suite apps are unreachable or the
/// suite cannot stay: App Translocation (Gatekeeper runs a quarantined app alone from a random
/// read-only path) or a read-only volume such as a mounted disk image.
enum InstallLocation {
    static let title = "Move KBioSuite to Applications"

    static func problem(bundlePath: String, volumeIsReadOnly: Bool?) -> String? {
        if bundlePath.contains("/AppTranslocation/") {
            return "macOS is running KBioSuite from a temporary, isolated location, so it cannot find the other suite apps. "
                + "Quit, drag the whole KBioSuite folder into your Applications folder, then open KBioSuite from there."
        }
        if volumeIsReadOnly == true {
            return "KBioSuite is running from a read-only volume, such as the disk image it came on. "
                + "Quit, drag the whole KBioSuite folder into your Applications folder, then open KBioSuite from there."
        }
        return nil
    }

    static func problem(bundleURL: URL) -> String? {
        let readOnly = try? bundleURL.resourceValues(forKeys: [.volumeIsReadOnlyKey]).volumeIsReadOnly
        return problem(bundlePath: bundleURL.path, volumeIsReadOnly: readOnly)
    }
}
