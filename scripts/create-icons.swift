import AppKit
import Foundation

let destination = URL(fileURLWithPath: CommandLine.arguments[1])
for (name, letters, color) in [("KBioSuite", "KB", NSColor(calibratedRed:0.18,green:0.34,blue:0.29,alpha:1)),
                               ("Kblast", "Kb", NSColor(calibratedRed:0.19,green:0.46,blue:0.37,alpha:1)),
                               ("Kannotate", "Ka", NSColor(calibratedRed:0.22,green:0.39,blue:0.54,alpha:1)),
                               ("Kmatch", "Km", NSColor(calibratedRed:0.43,green:0.28,blue:0.52,alpha:1))] {
    let directory = destination.appendingPathComponent(name + ".iconset")
    try FileManager.default.createDirectory(at: directory, withIntermediateDirectories:true)
    for size in [16,32,128,256,512] {
        for scale in [1,2] {
            let pixels=size*scale
            let image=NSImage(size:NSSize(width:pixels,height:pixels))
            image.lockFocus()
            color.setFill()
            NSBezierPath(roundedRect:NSRect(x:Double(pixels)*0.07,y:Double(pixels)*0.07,width:Double(pixels)*0.86,height:Double(pixels)*0.86),xRadius:Double(pixels)*0.19,yRadius:Double(pixels)*0.19).fill()
            let text=NSAttributedString(string:letters,attributes:[.font:NSFont.systemFont(ofSize:Double(pixels)*0.39,weight:.bold),.foregroundColor:NSColor.white])
            let extent=text.size();text.draw(at:NSPoint(x:(Double(pixels)-extent.width)/2,y:(Double(pixels)-extent.height)/2))
            image.unlockFocus()
            let bitmap=NSBitmapImageRep(data:image.tiffRepresentation!)!
            let filename="icon_\(size)x\(size)\(scale==2 ? "@2x" : "").png"
            try bitmap.representation(using:.png,properties:[:])!.write(to:directory.appendingPathComponent(filename))
        }
    }
    let process=Process();process.executableURL=URL(fileURLWithPath:"/usr/bin/iconutil")
    process.arguments=["-c","icns",directory.path,"-o",destination.appendingPathComponent(name+".icns").path]
    try process.run();process.waitUntilExit();guard process.terminationStatus==0 else { fatalError("iconutil failed") }
}
