import AppKit
import Foundation
let root=URL(fileURLWithPath:CommandLine.arguments[1])
let iconset=root.appendingPathComponent(".artifacts/Kstack.iconset")
try FileManager.default.createDirectory(at:iconset,withIntermediateDirectories:true)
for size in [16,32,128,256,512] {
    for factor in [1,2] {
        let px=size*factor,image=NSImage(size:NSSize(width:px,height:px))
        image.lockFocus()
        let scale=CGFloat(px)/1024
        let transform=NSAffineTransform();transform.scale(by:scale);transform.concat()
        NSColor(calibratedRed:0.15,green:0.36,blue:0.32,alpha:1).setFill()
        NSBezierPath(roundedRect:NSRect(x:70,y:70,width:884,height:884),xRadius:195,yRadius:195).fill()
        let colors=[NSColor(calibratedRed:0.91,green:0.95,blue:0.93,alpha:1),NSColor(calibratedRed:0.84,green:0.68,blue:0.46,alpha:1),NSColor(calibratedRed:0.60,green:0.71,blue:0.78,alpha:1)]
        for i in 0..<3 {
            let x:CGFloat=i==1 ? 294 : 230,y=CGFloat(584-i*192),path=NSBezierPath()
            path.move(to:NSPoint(x:x,y:y));path.line(to:NSPoint(x:x+390,y:y));path.line(to:NSPoint(x:x+500,y:y+80));path.line(to:NSPoint(x:x+390,y:y+160));path.line(to:NSPoint(x:x,y:y+160));path.close();colors[i].setFill();path.fill()
        }
        image.unlockFocus()
        let bitmap=NSBitmapImageRep(data:image.tiffRepresentation!)!
        let suffix=factor==2 ? "@2x" : ""
        try bitmap.representation(using:.png,properties:[:])!.write(to:iconset.appendingPathComponent("icon_\(size)x\(size)\(suffix).png"))
    }
}
