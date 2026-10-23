/* Excel 작업창. NCBI에는 직접 연결하지 않고 같은 기기의 blastid 서버만 호출한다. */
const DEFAULTS = {
  idHeader: "ID(Voucher)",
  seqHeader: "Seq",
  templateMode: "builtin",
  templateYaml: "",
  minIdent: 90,
  minCover: 80,
  altMargin: 2,
  batchSize: 10,
  email: "",
  offHours: false,
  backend: "auto",
  localDb: "ITS_RefSeq_Fungi",
  maxWait: 60,
};

const FILL = { strange: "#F8CBAD", ambiguous: "#FFF2CC", ok: "#E2EFDA" };

function norm(value) {
  return String(value ?? "").replace(/\s+/g, " ").trim().toLowerCase();
}

function parseRows(spec) {
  const rows = new Set();
  if (!spec || !spec.trim()) return null;
  for (const part of spec.split(",")) {
    const bit = part.trim();
    if (!bit) continue;
    if (bit.includes("-")) {
      const [a, b] = bit.split("-", 2).map((x) => parseInt(x, 10));
      if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
      for (let i = Math.min(a, b); i <= Math.max(a, b); i += 1) rows.add(i);
    } else {
      const n = parseInt(bit, 10);
      if (Number.isFinite(n)) rows.add(n);
    }
  }
  return rows;
}

function loadSettings() {
  let saved = null;
  try {
    const roaming = Office.context && Office.context.roamingSettings && Office.context.roamingSettings.get("blastid");
    if (roaming) saved = JSON.parse(roaming);
  } catch (_e) { /* 로밍 설정이 없으면 로컬로 */ }
  if (!saved) {
    try {
      const local = localStorage.getItem("blastid-settings");
      if (local) saved = JSON.parse(local);
    } catch (_e) { saved = null; }
  }
  return Object.assign({}, DEFAULTS, saved || {});
}

function saveSettings(settings) {
  const text = JSON.stringify(settings);
  try { localStorage.setItem("blastid-settings", text); } catch (_e) { /* 사생활 모드 */ }
  try {
    if (Office.context && Office.context.roamingSettings) {
      Office.context.roamingSettings.set("blastid", text);
      Office.context.roamingSettings.saveAsync(() => {});
    }
  } catch (_e) { /* 데스크톱이 아니면 로밍이 없을 수 있다 */ }
}

function readForm() {
  return {
    idHeader: document.getElementById("id-header").value.trim() || DEFAULTS.idHeader,
    seqHeader: document.getElementById("seq-header").value.trim() || DEFAULTS.seqHeader,
    templateMode: document.getElementById("template-mode").value,
    templateYaml: document.getElementById("template-file").dataset.yaml || "",
    minIdent: Number(document.getElementById("min-ident").value),
    minCover: Number(document.getElementById("min-cover").value),
    altMargin: Number(document.getElementById("alt-margin").value),
    batchSize: Number(document.getElementById("batch-size").value),
    email: document.getElementById("email").value.trim(),
    offHours: document.getElementById("off-hours").checked,
    backend: document.getElementById("backend").value,
    localDb: document.getElementById("local-db").value.trim() || "ITS_RefSeq_Fungi",
    maxWait: Number(document.getElementById("max-wait").value) || 60,
  };
}

function fillForm(settings) {
  document.getElementById("id-header").value = settings.idHeader;
  document.getElementById("seq-header").value = settings.seqHeader;
  document.getElementById("template-mode").value = settings.templateMode || "builtin";
  document.getElementById("min-ident").value = settings.minIdent;
  document.getElementById("min-cover").value = settings.minCover;
  document.getElementById("alt-margin").value = settings.altMargin;
  document.getElementById("batch-size").value = settings.batchSize;
  document.getElementById("email").value = settings.email || "";
  document.getElementById("off-hours").checked = Boolean(settings.offHours);
  document.getElementById("backend").value = settings.backend || "auto";
  document.getElementById("local-db").value = settings.localDb || "ITS_RefSeq_Fungi";
  document.getElementById("max-wait").value = settings.maxWait || 60;
  document.getElementById("template-file").dataset.yaml = settings.templateYaml || "";
  toggleYaml();
}

function toggleYaml() {
  const on = document.getElementById("template-mode").value === "yaml";
  document.getElementById("yaml-label").hidden = !on;
}

function columnIndex(headers, wanted) {
  const names = wanted.map(norm).filter(Boolean);
  for (let i = 0; i < headers.length; i += 1) {
    if (names.includes(norm(headers[i]))) return i;
  }
  return -1;
}

async function collectRows(settings) {
  const scope = document.querySelector("input[name=scope]:checked").value;
  const rangeSet = parseRows(document.getElementById("row-range").value);
  const boldOnly = document.getElementById("bold-only").checked;
  const onlyEmpty = document.getElementById("only-empty").checked;
  return Excel.run(async (context) => {
    const sheet = context.workbook.worksheets.getActiveWorksheet();
    const used = sheet.getUsedRange();
    used.load("rowCount, columnCount");
    let selection = null;
    if (scope === "selection") {
      selection = context.workbook.getSelectedRange();
      selection.load("rowIndex, rowCount");
    }
    await context.sync();
    const rows = Math.max(used.rowCount, 2);
    const cols = Math.max(used.columnCount, 1);
    const block = sheet.getRangeByIndexes(0, 0, rows, cols);
    block.load("values");
    await context.sync();
    const values = block.values;
    const headers = (values[0] || []).map((v) => String(v ?? "").trim());
    const idCol = columnIndex(headers, [settings.idHeader, "ID(Voucher)", "ID", "Voucher"]);
    const seqCol = columnIndex(headers, [settings.seqHeader, "Seq", "Sequence", "서열"]);
    const nameCol = columnIndex(headers, ["Scientific Name", "학명"]);
    if (seqCol < 0) throw new Error("서열 머리글을 찾지 못했습니다. 설정의 서열 머리글을 확인하세요.");
    let wanted = [];
    for (let r = 1; r < values.length; r += 1) wanted.push(r + 1);
    if (scope === "selection" && selection) {
      const start = selection.rowIndex + 1;
      const end = start + selection.rowCount - 1;
      wanted = wanted.filter((row) => row >= start && row <= end && row > 1);
    }
    if (rangeSet) wanted = wanted.filter((row) => rangeSet.has(row));
    const boldFlags = new Map();
    if (boldOnly && idCol >= 0) {
      const cells = wanted.map((row) => {
        const cell = sheet.getCell(row - 1, idCol);
        cell.load("format/font/bold");
        return cell;
      });
      await context.sync();
      wanted.forEach((row, i) => boldFlags.set(row, cells[i].format.font.bold === true));
    }
    const out = [];
    for (const row of wanted) {
      const line = values[row - 1] || [];
      const bold = boldOnly ? boldFlags.get(row) === true : false;
      if (boldOnly && !bold) continue;
      const name = nameCol >= 0 ? line[nameCol] : "";
      if (onlyEmpty && String(name ?? "").trim()) continue;
      out.push({
        row,
        id: idCol >= 0 ? line[idCol] : `row${row}`,
        seq: line[seqCol],
        name: name ?? "",
        bold,
      });
    }
    return { headers, rows: out };
  });
}

function findOrAddColumn(sheet, headers, spec) {
  const names = [spec.header].concat(spec.aliases || []);
  const existing = columnIndex(headers, names);
  if (existing >= 0) return existing;
  if (spec.role === "input" || !spec.field) return -1;
  const col = headers.length;
  headers.push(spec.header);
  const cell = sheet.getCell(0, col);
  cell.values = [[spec.header]];
  cell.format.fill.color = "#1F4E78";
  cell.format.font.color = "white";
  cell.format.font.bold = true;
  return col;
}

async function writeResults(result) {
  const columns = result.columns || [];
  const updates = result.updates || [];
  if (!updates.length) return;
  await Excel.run(async (context) => {
    const sheet = context.workbook.worksheets.getActiveWorksheet();
    const width = Math.max(columns.length + 8, 21);
    const header = sheet.getRangeByIndexes(0, 0, 1, width);
    header.load("values");
    await context.sync();
    const headers = header.values[0].map((v) => String(v ?? "").trim());
    while (headers.length && !headers[headers.length - 1]) headers.pop();
    const index = {};
    for (const spec of columns) {
      if (!spec.field) continue;
      index[spec.field] = findOrAddColumn(sheet, headers, spec);
    }
    const resultCols = columns.filter((c) => c.role === "result" && index[c.field] >= 0).map((c) => index[c.field]);
    for (const update of updates) {
      const r = update.row - 1;
      const local = update.backend === "local" || String(update.rid || "").startsWith("local");
      const fallback = update.level === "strange" ? FILL.strange : update.level === "ambiguous" ? FILL.ambiguous : local ? "#D6E3F0" : FILL.ok;
      const fill = update.fill ? `#${String(update.fill).replace(/^#/, "")}` : fallback;
      for (const [field, value] of Object.entries(update.values || {})) {
        const col = index[field];
        if (col == null || col < 0 || field === "screenshot") continue;
        const cell = sheet.getCell(r, col);
        cell.values = [[value]];
        const fmt = update.formats && update.formats[field];
        if (fmt) cell.numberFormat = [[fmt]];
      }
      for (const col of resultCols) {
        sheet.getCell(r, col).format.fill.color = fill;
      }
      if (local && index.source >= 0) {
        sheet.getCell(r, index.source).format.fill.color = "#D6E3F0";
      }
      sheet.getRangeByIndexes(r, 0, 1, 1).format.rowHeight = 30;
    }
    sheet.freezePanes.freezeAt(sheet.getRange("D2"));
    await context.sync();
    let anchor = 3;
    for (const update of updates) {
      if (update.comment && index.comment >= 0) {
        const addr = sheet.getCell(update.row - 1, index.comment);
        addr.load("address");
        await context.sync();
        try {
          context.workbook.comments.add(addr.address, update.comment);
          await context.sync();
        } catch (_e) { /* 이미 메모가 있으면 값만 남긴다 */ }
      }
      if (!update.screenshot_png) continue;
      let shot = context.workbook.worksheets.getItemOrNullObject("스크린샷");
      await context.sync();
      if (shot.isNullObject) {
        shot = context.workbook.worksheets.add("스크린샷");
        shot.getCell(0, 0).values = [["BLAST 상위 15개"]];
        shot.getRange("A:A").format.columnWidth = 120;
      }
      const cap = shot.getCell(anchor - 1, 0);
      cap.values = [[`${update.row}행 RID ${update.rid || ""}`]];
      try {
        const shape = shot.shapes.addImage(update.screenshot_png);
        shape.left = 8;
        shape.top = anchor * 18;
      } catch (_e) { /* 옛 Excel은 그림 API가 없을 수 있다 */ }
      if (index.screenshot >= 0) {
        const link = sheet.getCell(update.row - 1, index.screenshot);
        link.values = [[`스크린샷 보기 (${update.row}행)`]];
        link.format.font.color = "#0563C1";
        link.format.font.underline = "Single";
        try {
          link.hyperlink = {
            textToDisplay: `스크린샷 보기 (${update.row}행)`,
            documentReference: `'스크린샷'!A${anchor}`,
          };
        } catch (_e) { /* 내부 링크를 지원하지 않는 빌드 */ }
      }
      anchor += 20;
      await context.sync();
    }
  });
}

function appendLog(event) {
  const list = document.getElementById("log");
  const item = document.createElement("li");
  item.textContent = event.message || "";
  if (event.url) {
    const link = document.createElement("a");
    link.href = event.url;
    link.target = "_blank";
    link.rel = "noreferrer";
    link.textContent = event.rid ? ` RID ${event.rid}` : " NCBI 결과";
    item.appendChild(link);
  }
  list.appendChild(item);
  list.scrollTop = list.scrollHeight;
}

async function pollJob(id) {
  let seen = 0;
  for (;;) {
    const res = await fetch(`/api/jobs/${id}`);
    if (!res.ok) throw new Error(`서버 응답 ${res.status}`);
    const body = await res.json();
    for (const event of (body.events || []).slice(seen)) appendLog(event);
    seen = (body.events || []).length;
    document.getElementById("status").textContent = body.status === "waiting"
      ? "한가한 시간까지 기다리는 중"
      : body.status === "running" ? "NCBI 검색 중" : body.status === "done" ? "완료" : body.status === "error" ? "오류" : body.status;
    if (body.status === "done") return body.result || {};
    if (body.status === "error") throw new Error(body.error || "작업 실패");
    await new Promise((resolve) => setTimeout(resolve, 1500));
  }
}

async function runBlast() {
  const button = document.getElementById("run");
  button.disabled = true;
  document.getElementById("log").replaceChildren();
  const settings = readForm();
  saveSettings(settings);
  try {
    const collected = await collectRows(settings);
    appendLog({ message: `보낼 행 ${collected.rows.length}개` });
    if (!collected.rows.length) {
      document.getElementById("status").textContent = "조건에 맞는 행이 없습니다.";
      return;
    }
    const res = await fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        rows: collected.rows,
        bold_only: false,
        only_empty: false,
        min_ident: settings.minIdent,
        min_cover: settings.minCover,
        alt_margin: settings.altMargin,
        batch_size: settings.batchSize,
        email: settings.email,
        off_hours: settings.offHours,
        backend: settings.backend,
        local_db: settings.localDb,
        max_wait: settings.maxWait,
        template_yaml: settings.templateMode === "yaml" ? settings.templateYaml : null,
      }),
    });
    if (!res.ok) throw new Error(`서버 응답 ${res.status}`);
    const started = await res.json();
    const result = await pollJob(started.id);
    await writeResults(result);
    appendLog({ message: `시트에 ${((result && result.updates) || []).length}행을 기록했습니다.` });
  } catch (err) {
    const message = err && err.message ? err.message : String(err);
    const seen = [...document.querySelectorAll("#log li")].some((item) => (item.textContent || "").includes(message));
    if (!seen) appendLog({ message });
    document.getElementById("status").textContent = "실패";
  } finally {
    button.disabled = false;
  }
}

async function showDisclaimer() {
  try {
    const res = await fetch("/api/disclaimer");
    const body = await res.json();
    document.getElementById("disclaimer-ko").textContent = body.ko || "";
    document.getElementById("disclaimer-en").textContent = body.en || "";
    document.getElementById("window").textContent = body.window || "";
    const usage = await fetch("/api/usage");
    const counted = await usage.json();
    document.getElementById("usage").textContent =
      `최근 24시간 원격 검색 ${counted.count} / ${counted.limit}. 한도를 넘기면 나머지는 로컬 BLAST+입니다.`;
  } catch (_e) {
    document.getElementById("disclaimer-ko").textContent = "고지문을 서버에서 읽지 못했습니다. blastid serve 가 켜져 있는지 확인하세요.";
  }
}

function boot() {
  const settings = loadSettings();
  fillForm(settings);
  document.getElementById("template-mode").addEventListener("change", toggleYaml);
  document.getElementById("template-file").addEventListener("change", (event) => {
    const file = event.target.files && event.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      event.target.dataset.yaml = String(reader.result || "");
      const current = readForm();
      current.templateYaml = event.target.dataset.yaml;
      current.templateMode = "yaml";
      saveSettings(current);
    };
    reader.readAsText(file, "utf-8");
  });
  document.getElementById("save-settings").addEventListener("click", () => {
    saveSettings(readForm());
    document.getElementById("status").textContent = "설정을 저장했습니다.";
  });
  document.getElementById("run").addEventListener("click", runBlast);
  showDisclaimer();
}

Office.onReady(() => {
  if (typeof Excel === "undefined") {
    document.getElementById("status").textContent = "Excel 안에서 이 창을 여세요. 브라우저에서는 인증서 신뢰 확인용입니다.";
  }
  boot();
});
