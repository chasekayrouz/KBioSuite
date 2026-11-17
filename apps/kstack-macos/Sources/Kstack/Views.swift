import SwiftUI
import UniformTypeIdentifiers

private let forest = Color(red: 0.15, green: 0.37, blue: 0.32)

struct ContentView: View {
    @ObservedObject var model: AppModel
    @State private var showAnalysis = false
    var body: some View {
        ZStack {
            Group {
                if model.document.project != nil { WebEditor(model: model) }
                else { NewStackView(model: model) }
            }
            .disabled(model.busy || model.smokeRunning).allowsHitTesting(!model.busy && !model.smokeRunning)
            if model.busy || model.smokeRunning {
                Color.black.opacity(0.07).ignoresSafeArea()
                VStack(spacing: 14) {
                    ProgressView().controlSize(.large)
                    Text(model.status).font(.system(size: 13, weight: .medium))
                    Text("Your source files remain unchanged.").font(.caption).foregroundStyle(.secondary)
                }.padding(30).background(.regularMaterial, in: RoundedRectangle(cornerRadius: 14)).shadow(radius: 16)
            }
        }
        .frame(minWidth: 1080, minHeight: 720)
        .toolbar {
            ToolbarItemGroup(placement: .navigation) {
                Button { Task { await model.newStack() } } label: { Label("New Stack", systemImage: "plus") }.disabled(model.busy || model.smokeRunning)
                Button { model.browseOpen() } label: { Label("Open", systemImage: "folder") }.disabled(model.busy || model.smokeRunning)
            }
            ToolbarItemGroup(placement: .primaryAction) {
                if model.document.project != nil {
                    Button { showAnalysis = true } label: { Label("Homology", systemImage: "point.3.connected.trianglepath.dotted") }.disabled(model.busy || model.smokeRunning)
                    Button { Task { await model.save() } } label: { Label("Save", systemImage: "square.and.arrow.down") }.disabled(model.busy || model.smokeRunning)
                }
            }
        }
        .sheet(isPresented: $showAnalysis) {
            VStack(alignment: .leading, spacing: 19) {
                Text("Recompute homologs").font(.title2.weight(.semibold))
                Text("Search the stored proteins at new thresholds. Your layout and individual gene colors are preserved.").font(.callout).foregroundStyle(.secondary)
                HStack {
                    PercentField(label: "Amino-acid identity", value: $model.recomputeIdentity)
                    PercentField(label: "Shorter-protein coverage", value: $model.recomputeCoverage)
                }
                Text("Homolog groups use direct passing relationships and their transitive connections. Connector display filtering is controlled separately in the figure.").font(.caption).foregroundStyle(.secondary)
                HStack { Spacer(); Button("Cancel") { showAnalysis = false }; Button("Recompute") { showAnalysis = false; Task { await model.recompute() } }.buttonStyle(.borderedProminent).tint(forest) }
            }.padding(28).frame(width: 520)
        }
        .alert("Kstack", isPresented: Binding(get: { model.errorMessage != nil }, set: { if !$0 { model.errorMessage = nil } })) {
            Button("OK") { model.errorMessage = nil }
        } message: { Text(model.errorMessage ?? "") }
        .background(WindowHook(model: model))
    }
}

struct PercentField: View {
    let label: String
    @Binding var value: Double
    var body: some View {
        VStack(alignment: .leading, spacing: 6) {
            Text(label).font(.caption).foregroundStyle(.secondary)
            HStack(spacing: 4) {
                TextField(label, value: $value, format: .number.precision(.fractionLength(0...1))).textFieldStyle(.roundedBorder).frame(width: 70)
                Text("%").foregroundStyle(.secondary)
            }
        }
    }
}

struct NewStackView: View {
    @ObservedObject var model: AppModel
    @State private var targeted = false
    var body: some View {
        HStack(spacing: 0) {
            VStack(alignment: .leading, spacing: 22) {
                HStack(spacing: 10) {
                    Image(systemName: "arrow.right.square.fill").font(.system(size: 27)).foregroundStyle(forest)
                    Text("Kstack").font(.system(size: 23, weight: .semibold))
                }
                Text("A clear view of\nrelated gene clusters.").font(.system(size: 22, weight: .medium)).lineSpacing(5)
                Text("Import GenBanks, compare their proteins, and build an editable figure you can take anywhere.").font(.callout).foregroundStyle(.secondary).lineSpacing(4)
                Divider().padding(.vertical, 4)
                HStack { Text("RECENT STACKS").font(.system(size: 10, weight: .bold)).tracking(1.5).foregroundStyle(.secondary); Spacer() }
                ScrollView {
                    VStack(alignment: .leading, spacing: 12) {
                        if model.recent.isEmpty { Text("Saved stacks will appear here.").font(.caption).foregroundStyle(.tertiary) }
                        ForEach(model.recent) { entry in
                            HStack(alignment: .top) {
                                Button { Task { await model.open(entry.resolvedURL()) } } label: {
                                    HStack(alignment: .top, spacing: 9) {
                                        Image(systemName: entry.exists ? "doc.richtext" : "exclamationmark.triangle").foregroundStyle(entry.exists ? forest : .orange)
                                        VStack(alignment: .leading, spacing: 4) {
                                            Text(entry.name).font(.system(size: 12, weight: .medium)).lineLimit(2)
                                            Text(entry.exists ? entry.resolvedURL().deletingLastPathComponent().lastPathComponent : "File missing — locate with Open").font(.system(size: 10)).foregroundStyle(.secondary).lineLimit(1)
                                        }
                                    }.frame(maxWidth: .infinity, alignment: .leading)
                                }.buttonStyle(.plain).disabled(!entry.exists || model.busy).help(entry.path)
                                Button { model.removeRecent(entry) } label: { Image(systemName: "xmark").font(.system(size: 9)) }.buttonStyle(.plain).foregroundStyle(.tertiary).help("Remove recent reference")
                            }
                        }
                    }.padding(.vertical, 1)
                }
                Spacer(minLength: 0)
                Button("Open a saved HTML stack…") { model.browseOpen() }.buttonStyle(.link)
                Text("v\(Bundle.main.object(forInfoDictionaryKey: "CFBundleShortVersionString") as? String ?? "") · Local analysis · Portable HTML").font(.system(size: 10)).foregroundStyle(.tertiary)
            }.padding(28).frame(width: 280).background(Color(red: 0.94, green: 0.96, blue: 0.95))
            VStack(alignment: .leading, spacing: 18) {
                VStack(alignment: .leading, spacing: 8) {
                    Text("NEW STACK").font(.system(size: 10, weight: .bold)).tracking(1.8).foregroundStyle(forest)
                    Text("Start with your GenBank files").font(.system(size: 27, weight: .medium))
                    Text("1–50 files · One record per file · Up to 100 kb per record").font(.callout).foregroundStyle(.secondary)
                }
                VStack(spacing: 11) {
                    Image(systemName: "doc.badge.plus").font(.system(size: 30, weight: .light)).foregroundStyle(forest)
                    Text("Drop GenBank files here").font(.system(size: 15, weight: .medium))
                    Button("Browse / Add Files…") { model.addFiles() }.buttonStyle(.bordered).tint(forest)
                }.frame(maxWidth: .infinity).padding(23)
                    .background(targeted ? forest.opacity(0.12) : Color.white, in: RoundedRectangle(cornerRadius: 10))
                    .overlay(RoundedRectangle(cornerRadius: 10).strokeBorder(targeted ? forest : Color.gray.opacity(0.25), style: StrokeStyle(lineWidth: 1, dash: [5,4])))
                    .onDrop(of: [.fileURL], isTargeted: $targeted) { providers in
                        let group = DispatchGroup(); let collector = DropCollector()
                        for provider in providers {
                            group.enter(); provider.loadItem(forTypeIdentifier: UTType.fileURL.identifier, options: nil) { value, _ in
                                if let data = value as? Data, let url = URL(dataRepresentation: data, relativeTo: nil) { collector.add(url) }
                                else if let url = value as? URL { collector.add(url) }
                                group.leave()
                            }
                        }
                        group.notify(queue: .main) { Task { await model.importURLs(collector.urls) } }; return true
                    }
                ImportTable(model: model).frame(maxHeight: .infinity)
                VStack(alignment: .leading, spacing: 15) {
                    HStack(spacing: 28) {
                        PercentField(label: "Amino-acid identity", value: $model.identity)
                        PercentField(label: "Shorter-protein coverage", value: $model.coverage)
                        Spacer()
                        VStack(alignment: .leading, spacing: 7) {
                            Toggle("Preserve Input Order / Orientation", isOn: $model.preserve)
                            Toggle("Enable Similarity Connectors", isOn: $model.connectors)
                        }.font(.system(size: 11))
                    }
                    HStack {
                        VStack(alignment: .leading, spacing: 5) {
                            Text("OUTPUT HTML").font(.system(size: 9, weight: .bold)).tracking(1).foregroundStyle(.secondary)
                            HStack {
                                Text(model.outputURL?.path ?? "Choose a destination for your stack").font(.system(size: 11)).lineLimit(1).truncationMode(.middle).foregroundStyle(model.outputURL == nil ? .secondary : .primary)
                                Button("Choose…") { model.chooseOutput() }.controlSize(.small)
                            }
                            Text("Defaults to the first input’s folder. Choose… to change the name or location.").font(.system(size: 10)).foregroundStyle(.secondary)
                        }
                        Spacer()
                        Button { Task { await model.generate() } } label: { Label("Generate Stack", systemImage: "arrow.right").padding(.horizontal, 10).padding(.vertical, 5) }
                            .buttonStyle(.borderedProminent).tint(forest).disabled(!model.canGenerate)
                    }
                }.padding(18).background(Color.white, in: RoundedRectangle(cornerRadius: 10)).overlay(RoundedRectangle(cornerRadius: 10).stroke(.gray.opacity(0.15)))
            }.padding(30).background(Color(red: 0.975, green: 0.979, blue: 0.986))
        }.environment(\.colorScheme, .light)
    }
}
private final class DropCollector: @unchecked Sendable {
    private let lock = NSLock(); private var values: [URL] = []
    func add(_ url: URL) { lock.lock(); values.append(url); lock.unlock() }
    var urls: [URL] { lock.lock(); defer { lock.unlock() }; return values.sorted { $0.path < $1.path } }
}

struct ImportTable: View {
    @ObservedObject var model: AppModel
    var body: some View {
        VStack(spacing: 0) {
            HStack {
                Text("INPUT FILES (\(model.inputs.count))").font(.system(size: 10, weight: .bold)).tracking(1).foregroundStyle(.secondary)
                Spacer(); Text("Edit names before generation").font(.system(size: 10)).foregroundStyle(.tertiary)
            }.padding(.bottom, 10)
            if model.inputs.isEmpty {
                Text("Added files appear here for validation and ordering.").font(.callout).foregroundStyle(.tertiary).frame(maxWidth: .infinity, maxHeight: .infinity)
            } else {
                ScrollView {
                    LazyVStack(spacing: 0) {
                        ForEach($model.inputs) { $entry in
                            HStack(alignment: .top, spacing: 10) {
                                VStack(spacing: 5) {
                                    Button { model.moveInput(entry.id, by: -1) } label: { Image(systemName: "chevron.up") }
                                    Button { model.moveInput(entry.id, by: 1) } label: { Image(systemName: "chevron.down") }
                                }.font(.system(size: 9)).buttonStyle(.plain).foregroundStyle(.secondary)
                                VStack(alignment: .leading, spacing: 5) {
                                    TextField("Display name", text: $entry.name).textFieldStyle(.roundedBorder).font(.system(size: 12))
                                    Text(URL(fileURLWithPath: entry.path).lastPathComponent).font(.system(size: 10)).foregroundStyle(.secondary).lineLimit(1).help(entry.path)
                                }.frame(minWidth: 175, maxWidth: .infinity)
                                VStack(alignment: .trailing, spacing: 5) {
                                    Text(entry.cluster.map { String(format: "%0.1f kb", ($0["length"].number ?? 0)/1000) } ?? "—")
                                    Text("\(entry.cluster?["cds"].array.count ?? 0) CDS").foregroundStyle(.secondary)
                                }.font(.system(size: 11)).frame(width: 58)
                                VStack(alignment: .leading, spacing: 5) {
                                    if let error = entry.error { Label(error, systemImage: "xmark.circle.fill").foregroundStyle(.red) }
                                    else if entry.circular {
                                        if entry.circularChoice.isEmpty {
                                            Text("Circular record — choose:").foregroundStyle(.orange)
                                            HStack(spacing: 4) {
                                                Button("Linearize at Origin") { entry.circularChoice = "linearize" }
                                                Button("Omit") { entry.circularChoice = "omit" }
                                            }.controlSize(.mini)
                                        } else {
                                            Text(entry.omitted ? "Omitted" : "Linearized at existing origin").foregroundStyle(entry.omitted ? .secondary : forest)
                                            Button(entry.omitted ? "Linearize at Origin" : "Omit") { entry.circularChoice = entry.omitted ? "linearize" : "omit" }.controlSize(.mini)
                                        }
                                    } else { Label(entry.warningCount > 0 ? "\(entry.warningCount) warnings" : "Valid", systemImage: entry.warningCount > 0 ? "exclamationmark.triangle" : "checkmark.circle.fill").foregroundStyle(entry.warningCount > 0 ? .orange : forest) }
                                    if entry.warningCount > 0 {
                                        Text(warnings(entry)).lineLimit(2).foregroundStyle(.secondary).help(warnings(entry))
                                    }
                                }.font(.system(size: 10)).frame(width: 215, alignment: .leading)
                                Button { model.inputs.removeAll { $0.id == entry.id } } label: { Image(systemName: "xmark") }.buttonStyle(.plain).font(.system(size: 10)).foregroundStyle(.tertiary).help("Remove file")
                            }.padding(.vertical, 11)
                            Divider()
                        }
                    }
                }
            }
        }
    }
    func warnings(_ entry: ImportEntry) -> String {
        guard let c = entry.cluster else { return "" }
        return (c["warnings"].array + c["cds"].array.flatMap { $0["warnings"].array }).compactMap(\.string).joined(separator: " ")
    }
}
