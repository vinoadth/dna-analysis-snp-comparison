(function () {
  const MAX_BARS = 16;

  function $(sel, root) {
    return (root || document).querySelector(sel);
  }

  function fmtPct(value) {
    if (value === null || value === undefined || Number.isNaN(Number(value))) return "—";
    return Number(value).toFixed(2) + "%";
  }

  function barRows(items, className) {
    const cleaned = items
      .map(function (item) {
        return {
          name: item.population || item.label || "unknown",
          percent: item.percent,
        };
      })
      .filter(function (item) {
        return item.percent !== null && item.percent !== undefined;
      })
      .slice(0, MAX_BARS);
    if (!cleaned.length) return '<p class="text-secondary mb-0">No estimates for this comparison.</p>';
    const max = Math.max.apply(
      null,
      cleaned.map(function (item) {
        return Math.abs(Number(item.percent)) || 0;
      }).concat([1])
    );
    return (
      '<div class="d-grid gap-2">' +
      cleaned
        .map(function (item) {
          const width = Math.max(0, Math.min(100, (Math.abs(Number(item.percent)) / max) * 100));
          return (
            '<div class="row align-items-center g-2">' +
            '<div class="col-4 col-md-3 text-truncate small" title="' +
            escapeHtml(item.name) +
            '">' +
            escapeHtml(item.name) +
            "</div>" +
            '<div class="col">' +
            '<div class="progress" role="progressbar" aria-valuenow="' +
            Math.round(width) +
            '" aria-valuemin="0" aria-valuemax="100" style="height: 10px">' +
            '<div class="progress-bar ' +
            className +
            '" style="width:' +
            width +
            '%"></div></div></div>' +
            '<div class="col-auto font-monospace small">' +
            fmtPct(item.percent) +
            "</div></div>"
          );
        })
        .join("") +
      "</div>"
    );
  }

  function escapeHtml(text) {
    return String(text)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function chromOrder(chrom) {
    var key = String(chrom).replace(/^chr/i, "").toUpperCase();
    if (key === "X") return 23;
    if (key === "Y") return 24;
    if (key === "MT" || key === "M") return 25;
    var n = parseInt(key, 10);
    return Number.isFinite(n) ? n : 100;
  }

  function chromBars(counts) {
    const entries = Object.keys(counts || {})
      .map(function (chrom) {
        return { chrom: chrom, n: counts[chrom] };
      })
      .sort(function (a, b) {
        return chromOrder(a.chrom) - chromOrder(b.chrom);
      });
    if (!entries.length) return '<p class="text-secondary mb-0">No chromosome counts.</p>';
    const max = Math.max.apply(
      null,
      entries.map(function (item) {
        return item.n;
      })
    ) || 1;
    return (
      '<div class="d-grid gap-2">' +
      entries
        .map(function (item) {
          const width = (item.n / max) * 100;
          return (
            '<div class="row align-items-center g-2">' +
            '<div class="col-4 col-md-3 text-truncate small">chr ' +
            escapeHtml(item.chrom) +
            "</div>" +
            '<div class="col">' +
            '<div class="progress" role="progressbar" aria-valuenow="' +
            Math.round(width) +
            '" aria-valuemin="0" aria-valuemax="100" style="height: 10px">' +
            '<div class="progress-bar" style="width:' +
            width +
            '%"></div></div></div>' +
            '<div class="col-auto font-monospace small">' +
            item.n.toLocaleString() +
            "</div></div>"
          );
        })
        .join("") +
      "</div>"
    );
  }

  function fmtHgCell(cell) {
    if (!cell || !cell.n_called) return "—";
    const n = Number(cell.n || 0);
    const denom = Number(cell.n_called);
    if (cell.percent === null || cell.percent === undefined) return n + "/" + denom;
    return n + "/" + denom + " (" + Number(cell.percent).toFixed(0) + "%)";
  }

  const HG_STATUS_LABEL = {
    derived: "derived (yes)",
    ancestral: "ancestral (no)",
    "no-call": "no-call (SNP missing)",
    conflict: "conflict (ignore)",
    het: "het (unclear)",
    mismatch: "mismatch (unexpected allele)",
  };

  function hgStatusBadge(status) {
    const label = HG_STATUS_LABEL[status] || status;
    const cls = {
      derived: "text-bg-success",
      ancestral: "text-bg-secondary",
      "no-call": "text-bg-light",
      conflict: "text-bg-danger",
      het: "text-bg-warning",
      mismatch: "text-bg-danger",
    }[status] || "text-bg-light";
    return '<span class="badge ' + cls + '">' + escapeHtml(label) + "</span>";
  }

  function hgGlossary() {
    return (
      '<dl class="row small mb-3 p-3 bg-body-secondary rounded border">' +
      '<dt class="col-sm-2 font-monospace">derived</dt><dd class="col-sm-10">yes — this file has the mutation that defines that haplogroup</dd>' +
      '<dt class="col-sm-2 font-monospace">ancestral</dt><dd class="col-sm-10">no — this file has the older allele, so that haplogroup is ruled out</dd>' +
      '<dt class="col-sm-2 font-monospace">no-call</dt><dd class="col-sm-10">that defining SNP is missing or unreadable in this VCF</dd>' +
      '<dt class="col-sm-2 font-monospace">conflict</dt><dd class="col-sm-10">markers disagree (child looks yes, parent is no) — do not treat as a call</dd>' +
      '<dt class="col-sm-2 font-monospace">GQ</dt><dd class="col-sm-10">genotype quality from the VCF; higher is more confident (scale differs by file)</dd>' +
      '<dt class="col-sm-2 font-monospace">DP</dt><dd class="col-sm-10">sequencing read depth at that site; usually missing on SNP-array files</dd>' +
      '<dt class="col-sm-2 font-monospace">IGC</dt><dd class="col-sm-10">Illumina GenCall (0–1) on array VCFs; low IGC is weaker support</dd>' +
      "</dl>"
    );
  }

  function fmtQc(value, digits) {
    if (value === null || value === undefined || value === "") return "—";
    if (typeof value === "number" && digits !== undefined) return value.toFixed(digits);
    return String(value);
  }

  function fmtHgRsid(row) {
    const rsid = String(row.rsid || "").trim();
    if (/^rs\d/i.test(rsid)) {
      return '<span class="font-monospace">' + escapeHtml(rsid) + "</span>";
    }
    return (
      '<span class="text-secondary" title="No dbSNP rsID in the catalog; scored by chromosome position">' +
      "—</span>"
    );
  }

  function fmtHgAlleles(row) {
    const anc = row.ancestral || "—";
    const der = row.derived || "—";
    return (
      '<span class="font-monospace">' +
      escapeHtml(anc) +
      " → " +
      escapeHtml(der) +
      "</span>"
    );
  }

  function fmtHgObserved(row) {
    const obs = row.observed;
    const gt = row.genotype;
    if (obs == null && (!gt || gt === ".")) return "—";
    let text = obs != null ? String(obs) : "";
    if (gt && gt !== ".") {
      text = (text ? text + " " : "") + "(GT " + gt + ")";
    }
    return '<span class="font-monospace">' + escapeHtml(text || "—") + "</span>";
  }

  function hgMarkerSectionId(label) {
    return "hg-markers-" + String(label).toLowerCase().replace(/[^a-z0-9]+/g, "-");
  }

  function hgMarkerSearchHay(row) {
    return [
      row.haplogroup,
      row.marker,
      row.rsid,
      row.chrom,
      row.pos,
      row.ancestral,
      row.derived,
      row.observed,
      row.genotype,
      row.status,
    ]
      .filter(function (part) {
        return part != null && part !== "";
      })
      .join(" ")
      .toLowerCase();
  }

  function markerDetailTable(markers, label) {
    if (!markers || !markers.length) return "";
    const sectionId = hgMarkerSectionId(label);
    const hasDp = markers.some(function (m) { return m.dp != null; });
    const hasIgc = markers.some(function (m) { return m.igc != null; });
    const hasGq = markers.some(function (m) { return m.gq != null; });
    const hasQual = markers.some(function (m) { return m.qual != null; });
    const head =
      "<th>Haplogroup</th><th>Marker</th><th>rsID</th><th>Position</th><th>Defining alleles</th><th>Observed</th><th>This sample</th>" +
      (hasGq ? "<th>GQ</th>" : "") +
      (hasDp ? "<th>DP</th>" : "") +
      (hasIgc ? "<th>IGC</th>" : "") +
      (hasQual ? "<th>QUAL</th>" : "");
    const filters = [
      { id: "all", label: "All" },
      { id: "derived", label: "Derived" },
      { id: "lineage", label: "Derived + no-call" },
      { id: "ancestral", label: "Ancestral" },
      { id: "no-call", label: "No-call" },
      { id: "issues", label: "Issues" },
    ];
    const filterPills = filters
      .map(function (item, i) {
        return (
          '<button type="button" class="btn btn-sm hg-marker-filter-btn' +
          (i === 0 ? " btn-primary" : " btn-outline-secondary") +
          '" data-section="' +
          escapeHtml(sectionId) +
          '" data-filter="' +
          escapeHtml(item.id) +
          '">' +
          escapeHtml(item.label) +
          "</button>"
        );
      })
      .join("");
    const body = markers
      .map(function (row) {
        const status = row.status || "no-call";
        const pos =
          (row.chrom || "—") +
          ":" +
          (row.pos != null ? Number(row.pos).toLocaleString() : "—");
        return (
          '<tr data-hg-status="' +
          escapeHtml(status) +
          '" data-hg-search="' +
          escapeHtml(hgMarkerSearchHay(row)) +
          '"><td>' +
          escapeHtml(row.haplogroup) +
          "</td><td class=\"font-monospace\">" +
          escapeHtml(row.marker || "—") +
          "</td><td>" +
          fmtHgRsid(row) +
          '</td><td class="font-monospace">' +
          escapeHtml(pos) +
          "</td><td>" +
          fmtHgAlleles(row) +
          "</td><td>" +
          fmtHgObserved(row) +
          "</td><td>" +
          hgStatusBadge(status) +
          "</td>" +
          (hasGq ? '<td class="num">' + fmtQc(row.gq) + "</td>" : "") +
          (hasDp ? '<td class="num">' + fmtQc(row.dp) + "</td>" : "") +
          (hasIgc ? '<td class="num">' + fmtQc(row.igc, 2) + "</td>" : "") +
          (hasQual ? '<td class="num">' + fmtQc(row.qual) + "</td>" : "") +
          "</tr>"
        );
      })
      .join("");
    return (
      '<section class="hg-marker-section mt-4" id="' +
      escapeHtml(sectionId) +
      '" data-hg-filter="all">' +
      '<h3 class="h6 text-secondary">' +
      escapeHtml(label) +
      " defining SNP markers</h3>" +
      '<p class="small text-secondary mb-2">Backbone markers checked in this VCF: rsID (when catalogued), physical position, expected ancestral → derived allele, and your call. Derived = yes for that haplogroup.</p>' +
      '<div class="d-flex flex-wrap gap-2 align-items-center mb-2">' +
      '<input class="form-control form-control-sm hg-marker-search" style="max-width:16rem" type="search" data-section="' +
      escapeHtml(sectionId) +
      '" placeholder="Search rsID, marker, haplogroup" />' +
      '<span class="small text-secondary hg-marker-count" data-section="' +
      escapeHtml(sectionId) +
      '"></span>' +
      "</div>" +
      '<div class="d-flex flex-wrap gap-1 mb-2 hg-marker-filters" data-section="' +
      escapeHtml(sectionId) +
      '">' +
      filterPills +
      "</div>" +
      '<div class="table-responsive snp-table-wrap"><table class="table table-sm table-striped table-hover align-middle hg-table hg-marker-table"><thead><tr>' +
      head +
      "</tr></thead><tbody>" +
      body +
      "</tbody></table></div></section>"
    );
  }

  function hgMarkerMatchesFilter(status, filterId) {
    if (filterId === "all") return true;
    if (filterId === "derived") return status === "derived";
    if (filterId === "lineage") return status === "derived" || status === "no-call";
    if (filterId === "ancestral") return status === "ancestral";
    if (filterId === "no-call") return status === "no-call";
    if (filterId === "issues") {
      return status === "conflict" || status === "het" || status === "mismatch";
    }
    return true;
  }

  function applyHgMarkerFilter(section) {
    const filterId = section.getAttribute("data-hg-filter") || "all";
    const search = section.querySelector(".hg-marker-search");
    const q = search && search.value ? search.value.trim().toLowerCase() : "";
    const rows = section.querySelectorAll("tbody tr");
    let shown = 0;
    rows.forEach(function (row) {
      const status = row.getAttribute("data-hg-status") || "no-call";
      const hay = row.getAttribute("data-hg-search") || "";
      const match =
        hgMarkerMatchesFilter(status, filterId) && (!q || hay.indexOf(q) !== -1);
      row.classList.toggle("d-none", !match);
      if (match) shown += 1;
    });
    const count = section.querySelector(".hg-marker-count");
    if (count) {
      count.textContent =
        shown.toLocaleString() +
        " of " +
        rows.length.toLocaleString() +
        " markers shown";
    }
  }

  function bindHgMarkerFilters(root) {
    if (!root) return;
    root.querySelectorAll(".hg-marker-section").forEach(function (section) {
      applyHgMarkerFilter(section);
      section.querySelectorAll(".hg-marker-filter-btn").forEach(function (btn) {
        btn.addEventListener("click", function () {
          section.setAttribute("data-hg-filter", btn.getAttribute("data-filter") || "all");
          section.querySelectorAll(".hg-marker-filter-btn").forEach(function (peer) {
            const active = peer === btn;
            peer.classList.toggle("btn-primary", active);
            peer.classList.toggle("btn-outline-secondary", !active);
          });
          applyHgMarkerFilter(section);
        });
      });
      const search = section.querySelector(".hg-marker-search");
      if (search) {
        search.addEventListener("input", function () {
          applyHgMarkerFilter(section);
        });
      }
    });
  }

  function ancestry5Card(block) {
    if (!block || !block.available) {
      return "";
    }
    const rows = (block.estimates || [])
      .map(function (row) {
        return (
          "<tr><td>" +
          escapeHtml(row.population || "") +
          '</td><td class="num">' +
          (row.percent == null ? "—" : Number(row.percent).toFixed(1) + "%") +
          '</td><td class="num">' +
          Number(row.n_snps || 0).toLocaleString() +
          '</td><td class="num">' +
          (row.mean_ibs == null ? "—" : Number(row.mean_ibs).toFixed(3)) +
          "</td></tr>"
        );
      })
      .join("");
    return card(
      "Deep ancestry (qpAdm-style 5-source)",
      '<p class="small text-secondary mb-2">Same overlapping HO SNPs as the 3-source bars, plus Anatolia_N (Turkey_N) and East_Asian (Dai). A second model that sums to 100% — not a caste call.</p>' +
        '<div class="table-responsive mb-2"><table class="table table-sm table-striped table-hover align-middle hg-table"><thead><tr>' +
        "<th>Source</th><th>Percentage</th><th>Overlapping SNPs</th><th>Mean IBS</th>" +
        "</tr></thead><tbody>" +
        rows +
        "</tbody></table></div>" +
        notesList(block.notes)
    );
  }

  function additionalCard(block) {
    const rows = ((block && block.estimates) || [])
      .map(function (row) {
        const cls = row.available ? "" : "table-warning";
        const evidence = String(row.evidence || "—")
          .split(/\s*;\s*/)
          .filter(Boolean)
          .map(escapeHtml)
          .join("<br>");
        return (
          "<tr class=\"" +
          cls +
          "\"><td>" +
          escapeHtml(row.topic || "") +
          "</td><td>" +
          escapeHtml(row.coverage_status || (row.available ? "partial" : "missing")) +
          "</td><td class=\"col-finding\">" +
          escapeHtml(row.finding || "") +
          '</td><td class="num">' +
          Number(row.n_snps || 0).toLocaleString() +
          " / " +
          Number(row.n_markers || 0).toLocaleString() +
          "</td><td class=\"small col-evidence\">" +
          (evidence || "—") +
          "</td><td class=\"small col-source\">" +
          escapeHtml(row.source || "") +
          "</td><td class=\"small col-note\">" +
          escapeHtml(row.message || "—") +
          "</td></tr>"
        );
      })
      .join("");
    return card(
      "Additional details",
      '<p class="small text-secondary mb-2">Published-marker overlays (not a diagnosis): blood type, celiac/IBD/autoimmune tags, pharmacogenomics, metabolism, appearance, APOE, and more. Coverage shows whether this file has the needed SNPs (ready / partial / missing).</p>' +
        '<div class="table-responsive"><table class="table table-sm table-striped table-hover align-top additional-table"><thead><tr>' +
        "<th>Topic</th><th>Coverage</th><th>Finding</th><th>SNPs used / listed</th><th>Evidence</th><th>Source</th><th>Note</th>" +
        "</tr></thead><tbody>" +
        (rows || "<tr><td colspan=\"7\">No additional details yet.</td></tr>") +
        "</tbody></table></div>" +
        notesList(block && block.notes)
    );
  }

  function diseaseCoverageBadge(status) {
    const key = String(status || "missing").toLowerCase();
    const labels = {
      ready: "Ready",
      partial: "Partial",
      low: "Low overlap",
      missing: "Missing",
    };
    const classes = {
      ready: "text-bg-secondary",
      partial: "text-bg-light border text-dark",
      low: "text-bg-light border text-muted",
      missing: "text-bg-light border text-muted",
    };
    return (
      '<span class="badge disease-coverage-badge ' +
      (classes[key] || classes.missing) +
      '">' +
      escapeHtml(labels[key] || labels.missing) +
      "</span>"
    );
  }

  function diseaseRelativeLabel(row) {
    if (row.trait_category === "cognitive") {
      return row.percentile == null
        ? '<span class="text-secondary small">Research trait</span>'
        : '<span class="text-secondary small">' +
            Number(row.percentile).toFixed(0) +
            "th pct</span>";
    }
    const level = String(row.relative_level || "");
    if (level === "elevated") {
      return '<span class="badge disease-badge-elevated">Elevated</span>';
    }
    if (level === "typical") {
      return '<span class="text-secondary">Typical</span>';
    }
    if (level === "uncalibrated") {
      return '<span class="text-secondary small">Uncalibrated</span>';
    }
    return "—";
  }

  function diseaseCard(block) {
    const estimates = (block && block.estimates) || [];
    const elevatedCount = estimates.filter(function (row) {
      return row.relative_level === "elevated";
    }).length;
    const rows = estimates
      .map(function (row) {
        const cls = row.relative_level === "elevated" ? "disease-elevated" : "";
        return (
          "<tr class=\"" +
          cls +
          "\"><td>" +
          escapeHtml(row.trait || "") +
          "</td><td>" +
          diseaseCoverageBadge(row.coverage_status) +
          "</td><td>" +
          diseaseRelativeLabel(row) +
          '</td><td class="num">' +
          (row.percentile == null ? "—" : Number(row.percentile).toFixed(0)) +
          '</td><td class="num text-secondary small">' +
          (row.score == null ? "—" : Number(row.score).toFixed(3)) +
          '</td><td class="num">' +
          Number(row.n_snps || 0).toLocaleString() +
          " / " +
          Number(row.n_score || 0).toLocaleString() +
          '</td><td class="num">' +
          (row.coverage_pct == null ? "—" : Number(row.coverage_pct).toFixed(1) + "%") +
          "</td><td class=\"small\">" +
          escapeHtml(row.pgs_id || "") +
          "</td><td class=\"small\">" +
          escapeHtml(row.message || "—") +
          "</td></tr>"
        );
      })
      .join("");
    const summary =
      elevatedCount > 0
        ? '<p class="small mb-2"><span class="badge disease-badge-elevated">' +
          elevatedCount +
          " elevated</span> — at or above the 90th percentile in a European reference (≥50% SNP overlap).</p>"
        : '<p class="small text-secondary mb-2">No elevated disease scores in this run. Amber highlighting appears only for ≥90th-percentile results with ≥50% SNP overlap.</p>';
    return card(
      "Published polygenic scores (research only)",
      summary +
        '<p class="small text-secondary mb-2">Research polygenic scores from the PGS Catalog — not a diagnosis or lifetime risk. Coverage badges show data quality (not genetic risk). Cognitive traits are shown without disease-risk highlighting. There is no single “all cancer” score.</p>' +
        '<div class="table-responsive"><table class="table table-sm table-striped table-hover align-middle hg-table disease-table"><thead><tr>' +
        "<th>Trait</th><th>Coverage</th><th>Relative level</th><th>Percentile</th><th>Raw score</th><th>SNPs used / in score</th><th>Overlap</th><th>PGS</th><th>Note</th>" +
        "</tr></thead><tbody>" +
        (rows || "<tr><td colspan=\"9\">No disease scores.</td></tr>") +
        "</tbody></table></div>" +
        notesList(block && block.notes)
    );
  }

  function communityRefRowClass(row) {
    const usesIndus = Boolean(row.ref_indus);
    const usesEastAsian = Boolean(row.ref_east_asian);
    if (usesIndus && row.aasi_in_range && row.steppe_in_range && row.indus_in_range) {
      return "table-success";
    }
    if (usesEastAsian && row.aasi_in_range && row.steppe_in_range && row.east_asian_in_range) {
      return "table-success";
    }
    if (row.steppe_in_range && row.aasi_in_range) {
      return "table-success";
    }
    if (row.steppe_in_range || row.indus_in_range || row.east_asian_in_range) {
      return "table-info";
    }
    if (row.percent >= 50) {
      return "table-warning";
    }
    return "";
  }

  function communityRefPanelCard(title, rows, notes) {
    const sample = rows[0] || {};
    const usesIndus = rows.some(function (row) {
      return row.ref_indus;
    });
    const usesEastAsian = rows.some(function (row) {
      return row.ref_east_asian;
    });
    const body = rows
      .map(function (row) {
        return (
          "<tr class=\"" +
          communityRefRowClass(row) +
          "\"><td>" +
          escapeHtml(row.population || "") +
          '</td><td class="num">' +
          escapeHtml(row.ref_aasi || "—") +
          '</td><td class="num">' +
          escapeHtml(row.ref_steppe || "—") +
          (usesIndus
            ? '</td><td class="num">' + escapeHtml(row.ref_indus || "—")
            : "") +
          (usesEastAsian
            ? '</td><td class="num">' + escapeHtml(row.ref_east_asian || "—")
            : "") +
          '</td><td class="num">' +
          (row.sample_aasi == null ? "—" : Number(row.sample_aasi).toFixed(1) + "%") +
          '</td><td class="num">' +
          (row.sample_steppe == null ? "—" : Number(row.sample_steppe).toFixed(1) + "%") +
          (usesIndus
            ? '</td><td class="num">' +
              (row.sample_indus == null ? "—" : Number(row.sample_indus).toFixed(1) + "%")
            : "") +
          (usesEastAsian
            ? '</td><td class="num">' +
              (row.sample_east_asian == null ? "—" : Number(row.sample_east_asian).toFixed(1) + "%")
            : "") +
          "</td><td>" +
          escapeHtml(row.ref_y || "—") +
          '</td><td class="num">' +
          Number(row.percent).toFixed(0) +
          "%</td></tr>"
        );
      })
      .join("");
    const intro = usesIndus
      ? "Fit compares AASI_Onge, Steppe_MLBA, and Indus_Periphery to published ranges. Not a caste call."
      : usesEastAsian
        ? "Fit compares ASI→AASI_Onge, ANI→Steppe_MLBA, and East Asian→East_Asian (5-source). Not a caste call."
        : "Fit compares AASI_Onge and Steppe_MLBA to published ranges. Not a caste call.";
    const head =
      "<th>Community</th><th>" +
      escapeHtml(sample.ref_aasi_label || "Ref AASI") +
      "</th><th>" +
      escapeHtml(sample.ref_steppe_label || "Ref Steppe") +
      "</th>" +
      (usesIndus ? "<th>" + escapeHtml(sample.ref_indus_label || "Ref Indus") + "</th>" : "") +
      (usesEastAsian
        ? "<th>" + escapeHtml(sample.ref_east_asian_label || "Ref E Asian") + "</th>"
        : "") +
      "<th>" +
      escapeHtml(sample.sample_aasi_label || "This AASI") +
      "</th><th>" +
      escapeHtml(sample.sample_steppe_label || "This Steppe") +
      "</th>" +
      (usesIndus ? "<th>" + escapeHtml(sample.sample_indus_label || "This Indus") + "</th>" : "") +
      (usesEastAsian
        ? "<th>" + escapeHtml(sample.sample_east_asian_label || "This E Asian") + "</th>"
        : "") +
      "<th>Ref Y</th><th>Fit</th>";
    return card(
      title,
      '<p class="small text-secondary mb-2">' +
        intro +
        " Row color: all key ranges match / Steppe or Indus matches / fit ≥ 50%.</p>" +
        '<div class="table-responsive"><table class="table table-sm table-striped table-hover align-middle hg-table"><thead><tr>' +
        head +
        "</tr></thead><tbody>" +
        body +
        "</tbody></table></div>" +
        notesList(notes)
    );
  }

  function communityRefCard(block) {
    if (!block || block.hidden) {
      return "";
    }
    if (!block.available) {
      return card(
        "Community reference ranges",
        '<p class="text-secondary mb-0">Need ancestry results and a matching community reference table under data/references/caste/.</p>' +
          notesList(block && block.notes)
      );
    }
    const rows = block.estimates || [];
    const panels = {};
    rows.forEach(function (row) {
      const key = row.panel || "Community reference ranges";
      if (!panels[key]) {
        panels[key] = [];
      }
      panels[key].push(row);
    });
    return Object.keys(panels)
      .map(function (title) {
        const panelRows = panels[title];
        const panelNotes = (block.notes || []).filter(function (note) {
          return String(note).indexOf(title) === 0 || String(note).indexOf(title + ":") === 0;
        });
        return communityRefPanelCard(title, panelRows, panelNotes.length ? panelNotes : null);
      })
      .join("");
  }

  function relatednessCard(block) {
    if (!block || !block.available) {
      if (!block || !block.other_filename) {
        return "";
      }
      return card(
        "Relatedness vs second VCF",
        '<p class="text-secondary mb-0">The second VCF was attached, but relatedness could not be estimated.</p>' +
          notesList(block.notes)
      );
    }
    const kinship = block.kinship == null ? "—" : Number(block.kinship).toFixed(3);
    const ibs = block.mean_ibs == null ? "—" : Number(block.mean_ibs).toFixed(3);
    const n = Number(block.n_snps || 0).toLocaleString();
    const rel = block.reliability || "—";
    const het =
      (block.het_rate_query == null ? "—" : Number(block.het_rate_query).toFixed(2)) +
      " / " +
      (block.het_rate_other == null ? "—" : Number(block.het_rate_other).toFixed(2));
    const matchBits = [];
    if (block.n_matched_pos) matchBits.push(Number(block.n_matched_pos).toLocaleString() + " by position");
    if (block.n_matched_rsid) matchBits.push(Number(block.n_matched_rsid).toLocaleString() + " by rsID");
    const body =
      '<p class="mb-3">Compared to <strong>' +
      escapeHtml(block.other_sample_id || block.other_filename || "second file") +
      "</strong></p>" +
      '<p class="fs-5 mb-2">' +
      escapeHtml(block.relationship || "unknown") +
      "</p>" +
      '<div class="row g-3 mb-2">' +
      kpi(kinship, "KING kinship") +
      kpi(ibs, "mean IBS") +
      kpi(n, "overlapping autosomal SNPs") +
      kpi(
        Number(block.ibs0 || 0).toLocaleString() +
          " / " +
          Number(block.ibs1 || 0).toLocaleString() +
          " / " +
          Number(block.ibs2 || 0).toLocaleString(),
        "IBS0 / IBS1 / IBS2"
      ) +
      "</div>" +
      '<div class="row g-3 mb-2">' +
      kpi(escapeHtml(rel), "call reliability") +
      kpi(het, "het rate (query / other)") +
      kpi(
        matchBits.length ? matchBits.join(" · ") : "—",
        "how SNPs were matched"
      ) +
      kpi(
        Number(block.n_qc_dropped || 0).toLocaleString(),
        "low-quality sites dropped"
      ) +
      "</div>" +
      relatednessTables(block) +
      notesList(block.notes);
    return card("Relatedness vs second VCF", body);
  }

  function relatednessBand(kinship) {
    if (kinship == null) return "";
    const k = Number(kinship);
    if (k >= 0.354) return "twin";
    if (k >= 0.177) return "first";
    if (k >= 0.088) return "second";
    if (k >= 0.044) return "cousin";
    if (k >= 0.022) return "distant";
    return "unrelated";
  }

  function relatednessTables(block) {
    const n = Number(block.n_snps || 0);
    const fmt = function (value) {
      if (value == null || Number.isNaN(Number(value))) return "—";
      return Number(value).toFixed(1) + "%";
    };
    const rows = [
      ["IBS0 (opposite homozygotes)", block.ibs0, block.ibs0_pct],
      ["IBS1 (one allele shared)", block.ibs1, block.ibs1_pct],
      ["IBS2 (both alleles shared)", block.ibs2, block.ibs2_pct],
      ["Mean IBS (allele sharing)", "—", block.mean_ibs == null ? null : 100 * Number(block.mean_ibs)],
      ["Estimated DNA shared (2 × kinship)", "—", block.shared_pct],
    ];
    const share = rows
      .map(function (row) {
        return (
          "<tr><td>" +
          escapeHtml(row[0]) +
          '</td><td class="num">' +
          (row[1] === "—" ? "—" : Number(row[1] || 0).toLocaleString()) +
          '</td><td class="num">' +
          fmt(row[2]) +
          "</td></tr>"
        );
      })
      .join("");
    const band = relatednessBand(block.kinship);
    const refs = [
      ["twin", "Same person / identical twin", "≥ 0.35", "~100%", "~0%"],
      ["first", "Parent–child or full sibling", "~0.25", "~50%", "≈0% parent–child"],
      ["second", "Second-degree (half-sib, uncle, grandparent)", "~0.13", "~25%", "low"],
      ["cousin", "Third-degree (first cousin)", "~0.06", "~12.5%", "moderate"],
      ["distant", "Fourth-degree / distant", "~0.03", "~6%", "higher"],
      ["unrelated", "Unrelated or very distant", "~0", "~0%", "highest"],
    ];
    const refBody = refs
      .map(function (row) {
        const active = row[0] === band;
        return (
          "<tr" +
          (active ? ' class="table-info"' : "") +
          "><td>" +
          escapeHtml(row[1]) +
          (active ? ' <span class="badge text-bg-info">this pair</span>' : "") +
          "</td><td>" +
          escapeHtml(row[2]) +
          "</td><td>" +
          escapeHtml(row[3]) +
          "</td><td>" +
          escapeHtml(row[4]) +
          "</td></tr>"
        );
      })
      .join("");
    return (
      '<h3 class="h6 text-secondary mt-3">Sharing on ' +
      n.toLocaleString() +
      " overlapping SNPs</h3>" +
      '<div class="table-responsive mb-3"><table class="table table-sm table-striped align-middle hg-table"><thead><tr>' +
      "<th>Metric</th><th>Count</th><th>Percentage</th></tr></thead><tbody>" +
      share +
      "</tbody></table></div>" +
      '<h3 class="h6 text-secondary">Typical ranges (KING)</h3>' +
      '<div class="table-responsive"><table class="table table-sm table-striped align-middle hg-table"><thead><tr>' +
      "<th>Relationship</th><th>Kinship</th><th>DNA shared</th><th>IBS0</th></tr></thead><tbody>" +
      refBody +
      "</tbody></table></div>"
    );
  }

  function haploCard(hg) {
    const showY = !!hg.available;
    const showMt = !!hg.mt_available;
    if (!showY && !showMt) {
      return card("Haplogroups", notesList(hg.notes) + notesList(hg.mt_notes));
    }
    const title = showY
      ? "Y haplogroups (R1a1 / M17 and others)"
      : "mtDNA haplogroups (M, R, U, and others)";
    let body = "";
    if (showY) {
      body += haploTable(hg, "y");
      body += markerDetailTable(hg.markers || [], "Y");
    }
    if (showMt) {
      if (showY) body += '<h2 class="h5 mt-4">mtDNA haplogroups (M, R, U, and others)</h2>';
      body += haploTable(
        { sample_best: hg.mt_sample_best, rows: hg.mt_rows || [] },
        "mt"
      );
      body += markerDetailTable(hg.mt_markers || [], "mtDNA");
    }
    if (showY || showMt) body += hgGlossary();
    body += notesList(hg.notes);
    body += notesList(hg.mt_notes);
    return card(title, body);
  }

  function hgPctClass(cell) {
    if (!cell || cell.percent === null || cell.percent === undefined) return "";
    const pct = Number(cell.percent);
    if (!(pct > 0)) return "";
    if (pct <= 25) return "table-warning";
    if (pct <= 50) return "table-info";
    if (pct <= 75) return "table-primary";
    return "table-success";
  }

  function hgPctLegend() {
    return (
      '<p class="small text-secondary mb-2">' +
      "Highlight is the AADR share in that community: " +
      '<span class="badge text-bg-warning">≤25%</span> ' +
      '<span class="badge text-bg-info">≤50%</span> ' +
      '<span class="badge text-bg-primary">≤75%</span> ' +
      '<span class="badge text-bg-success">&gt;75%</span>' +
      "</p>"
    );
  }

  function haploTable(block, kind) {
    const label = kind === "mt" ? "mtDNA" : "Y";
    const rows = (block && block.rows) || [];
    const best = block && block.sample_best
      ? '<p class="mb-3">Deepest derived ' +
        label +
        " marker in this VCF: <strong>" +
        escapeHtml(block.sample_best) +
        "</strong></p>"
      : '<p class="mb-3">No derived backbone ' +
        label +
        " marker in this VCF (or calls conflict).</p>";
    if (!rows.length) {
      return best + '<p class="text-secondary mb-0">No AADR haplogroup counts (need the .anno file).</p>';
    }
    const groups = Object.keys(rows[0].groups || {});
    const head =
      "<th>Haplogroup</th><th>Marker</th><th>This sample</th>" +
      groups
        .map(function (name) {
          return "<th>" + escapeHtml(name) + "</th>";
        })
        .join("");
    const body = rows
      .map(function (row) {
        const status = row.sample_status || "no-call";
        return (
          "<tr><td>" +
          escapeHtml(row.haplogroup) +
          "</td><td>" +
          escapeHtml(row.marker || "—") +
          "</td><td>" +
          hgStatusBadge(status) +
          "</td>" +
          groups
            .map(function (name) {
              const cell = (row.groups || {})[name];
              const extra = hgPctClass(cell);
              return (
                '<td class="num' +
                (extra ? " " + extra : "") +
                '">' +
                fmtHgCell(cell) +
                "</td>"
              );
            })
            .join("") +
          "</tr>"
        );
      })
      .join("");
    return (
      best +
      hgPctLegend() +
      '<div class="table-responsive"><table class="table table-sm table-striped table-hover align-middle hg-table"><thead><tr>' +
      head +
      "</tr></thead><tbody>" +
      body +
      "</tbody></table></div>"
    );
  }

  function notesList(notes) {
    if (!notes || !notes.length) return "";
    return (
      '<ul class="list-unstyled small text-secondary mb-0 mt-3">' +
      notes
        .map(function (note) {
          return "<li>" + escapeHtml(note) + "</li>";
        })
        .join("") +
      "</ul>"
    );
  }

  function snpRowHtml(row) {
    return (
      "<tr><td>" +
      escapeHtml(row.chrom) +
      '</td><td class="num">' +
      escapeHtml(row.pos) +
      "</td><td>" +
      escapeHtml(row.rsid) +
      "</td><td>" +
      escapeHtml(row.ref) +
      "</td><td>" +
      escapeHtml(row.alt) +
      "</td><td>" +
      escapeHtml(row.genotype) +
      "</td></tr>"
    );
  }

  function snpTableHtml(rows) {
    if (!rows || !rows.length) {
      return '<p class="text-secondary mb-0">No SNPs match this chromosome or search.</p>';
    }
    return (
      '<div class="table-responsive snp-table-wrap"><table class="table table-sm table-striped table-hover align-middle hg-table"><thead><tr>' +
      "<th>chrom</th><th>pos</th><th>rsid</th><th>ref</th><th>alt</th><th>GT</th>" +
      "</tr></thead><tbody>" +
      rows.map(snpRowHtml).join("") +
      "</tbody></table></div>"
    );
  }

  function snpListCard(vcf, sourceFilename) {
    const counts = vcf.chrom_counts || {};
    const chroms = Object.keys(counts).sort(function (a, b) {
      return chromOrder(a) - chromOrder(b);
    });
    const total = Number(vcf.n_catalog || vcf.n_snps || (vcf.preview || []).length);
    const pills = ['<button type="button" class="btn btn-sm btn-primary snp-chrom-btn" data-chrom="">All</button>']
      .concat(
        chroms.map(function (chrom) {
          return (
            '<button type="button" class="btn btn-sm btn-outline-primary snp-chrom-btn" data-chrom="' +
            escapeHtml(chrom) +
            '">chr ' +
            escapeHtml(chrom) +
            " <span class=\"badge text-bg-light\">" +
            Number(counts[chrom] || 0).toLocaleString() +
            "</span></button>"
          );
        })
      )
      .join("");
    const catalog = vcf.snp_catalog_id || sourceFilename || "";
    const download = catalog
      ? '<a class="btn btn-sm btn-outline-secondary" id="snp-download" href="/api/snps.tsv?file=' +
        encodeURIComponent(sourceFilename || catalog) +
        '">Download TSV</a>'
      : "";
    return card(
      "SNP list (all chromosomes)",
      '<p class="small text-secondary mb-2">' +
        total.toLocaleString() +
        " SNPs across " +
        chroms.length +
        " chromosomes. Filter by chromosome or search rsID / position.</p>" +
        '<div class="d-flex flex-wrap gap-3 align-items-center mb-2">' +
        '<input class="form-control form-control-sm" style="max-width:16rem" id="snp-search" type="search" placeholder="Search rsID or position" />' +
        download +
        '<span class="small text-secondary" id="snp-page-label"></span>' +
        "</div>" +
        '<div class="d-flex flex-wrap gap-1 mb-2" id="snp-chrom-filters">' +
        pills +
        "</div>" +
        '<div id="snp-table">' +
        snpTableHtml(vcf.preview || []) +
        "</div>" +
        '<div class="d-flex gap-2 mt-2">' +
        '<button type="button" class="btn btn-sm btn-outline-secondary" id="snp-prev">Previous</button>' +
        '<button type="button" class="btn btn-sm btn-outline-secondary" id="snp-next">Next</button>' +
        "</div>"
    );
  }

  function bindSnpList(root, payload) {
    const vcf = (payload && payload.vcf) || {};
    const filename = (payload && payload.source_filename) || "";
    const table = root.querySelector("#snp-table");
    if (!table) return;
    const state = {
      chrom: "",
      q: "",
      offset: 0,
      limit: 200,
      total: Number(vcf.n_catalog || 0),
    };
    const label = root.querySelector("#snp-page-label");
    const search = root.querySelector("#snp-search");
    const download = root.querySelector("#snp-download");
    const useApi = !!vcf.snp_catalog_id && !window.ANALYSIS_PAYLOAD;

    function localRows() {
      const rows = (vcf.preview || []).filter(function (row) {
        if (state.chrom && String(row.chrom) !== state.chrom) return false;
        if (!state.q) return true;
        const hay = [row.chrom, row.pos, row.rsid, row.ref, row.alt, row.genotype].join(" ").toLowerCase();
        return hay.indexOf(state.q) !== -1;
      });
      state.total = rows.length;
      return rows.slice(state.offset, state.offset + state.limit);
    }

    function syncDownload() {
      if (!download) return;
      const params = new URLSearchParams({ file: filename });
      if (state.chrom) params.set("chrom", state.chrom);
      download.href = "/api/snps.tsv?" + params.toString();
    }

    function paint(rows) {
      table.innerHTML = snpTableHtml(rows);
      const start = state.total ? state.offset + 1 : 0;
      const end = Math.min(state.offset + state.limit, state.total);
      if (label) {
        label.textContent = state.total
          ? "Showing " + start.toLocaleString() + "–" + end.toLocaleString() + " of " + state.total.toLocaleString()
          : "No SNPs";
      }
      const prev = root.querySelector("#snp-prev");
      const next = root.querySelector("#snp-next");
      if (prev) prev.disabled = state.offset <= 0;
      if (next) next.disabled = state.offset + state.limit >= state.total;
    }

    function refresh() {
      syncDownload();
      if (!useApi) {
        paint(localRows());
        return;
      }
      syncDownload();
      const params = new URLSearchParams({
        file: filename,
        chrom: state.chrom,
        q: state.q,
        offset: String(state.offset),
        limit: String(state.limit),
      });
      fetch("/api/snps?" + params.toString())
        .then(function (res) {
          return res.json();
        })
        .then(function (data) {
          state.total = Number(data.total || 0);
          paint(data.rows || []);
        })
        .catch(function () {
          paint(localRows());
        });
    }

    root.querySelectorAll(".snp-chrom-btn").forEach(function (btn) {
      btn.addEventListener("click", function () {
        state.chrom = btn.getAttribute("data-chrom") || "";
        state.offset = 0;
        root.querySelectorAll(".snp-chrom-btn").forEach(function (other) {
          other.className =
            "btn btn-sm " +
            (other === btn ? "btn-primary" : "btn-outline-primary") +
            " snp-chrom-btn";
        });
        refresh();
      });
    });
    if (search) {
      search.addEventListener("input", function () {
        state.q = search.value.trim().toLowerCase();
        state.offset = 0;
        refresh();
      });
    }
    const prev = root.querySelector("#snp-prev");
    const next = root.querySelector("#snp-next");
    if (prev) {
      prev.addEventListener("click", function () {
        state.offset = Math.max(0, state.offset - state.limit);
        refresh();
      });
    }
    if (next) {
      next.addEventListener("click", function () {
        if (state.offset + state.limit < state.total) {
          state.offset += state.limit;
          refresh();
        }
      });
    }
    refresh();
  }

  function card(title, body) {
    return (
      '<section class="card shadow-sm mb-3">' +
      '<div class="card-header fw-semibold">' +
      escapeHtml(title) +
      "</div>" +
      '<div class="card-body">' +
      body +
      "</div></section>"
    );
  }

  function kpi(value, label, valueClass) {
    return (
      '<div class="col-6 col-md-3">' +
      '<div class="card shadow-sm h-100"><div class="card-body py-3">' +
      '<div class="fs-4 fw-semibold ' +
      (valueClass || "") +
      '">' +
      value +
      '</div><div class="text-secondary small">' +
      label +
      "</div></div></div></div>"
    );
  }

  function resultTabs(payload) {
    const vcf = payload.vcf || {};
    const related = relatednessCard(payload.relatedness || {});
    const tabs = [
      {
        id: "ancestry",
        label: "Ancestry",
        body:
          card(
            "Deep ancestry (qpAdm-style 3-source)",
            barRows((payload.ancestry && payload.ancestry.estimates) || [], "ancestry") +
              notesList(payload.ancestry && payload.ancestry.notes)
          ) +
          ancestry5Card(payload.ancestry_5 || {}),
      },
      {
        id: "communities",
        label: "Communities",
        body:
          '<div class="row g-3">' +
          '<div class="col-lg-6">' +
          card(
            "Population mixture weights",
            barRows((payload.populations && payload.populations.estimates) || [], "") +
              notesList(payload.populations && payload.populations.notes)
          ) +
          "</div><div class=\"col-lg-6\">" +
          card(
            "Caste / community weights",
            barRows((payload.caste && payload.caste.estimates) || [], "") +
              notesList(payload.caste && payload.caste.notes)
          ) +
          "</div></div>" +
          communityRefCard(payload.community_ref || {}),
      },
      {
        id: "haplogroups",
        label: "Haplogroups",
        body: haploCard(payload.haplogroups || {}),
      },
    ];
    tabs.push({
      id: "disease",
      label: "Disease (research)",
      body: diseaseCard(payload.disease || {}),
    });
    tabs.push({
      id: "additional",
      label: "Additional details",
      body: additionalCard(payload.additional || {}),
    });
    if (related) {
      tabs.push({ id: "relatedness", label: "Relatedness", body: related });
    }
    tabs.push({
      id: "hominin",
      label: "Hominin",
      body: card(
        "Hominin comparison",
        barRows((payload.hominin && payload.hominin.estimates) || [], "hominin") +
          notesList(payload.hominin && payload.hominin.notes)
      ),
    });
    tabs.push({
      id: "snps",
      label: "SNPs",
      body:
        card("SNPs by chromosome", chromBars(vcf.chrom_counts || {})) +
        snpListCard(vcf, payload.source_filename || ""),
    });
    const nav = tabs
      .map(function (tab, i) {
        return (
          '<li class="nav-item" role="presentation">' +
          '<button class="nav-link' +
          (i === 0 ? " active" : "") +
          '" id="tab-' +
          tab.id +
          '-btn" data-bs-toggle="tab" data-bs-target="#tab-' +
          tab.id +
          '" type="button" role="tab" aria-controls="tab-' +
          tab.id +
          '" aria-selected="' +
          (i === 0 ? "true" : "false") +
          '">' +
          escapeHtml(tab.label) +
          "</button></li>"
        );
      })
      .join("");
    const panes = tabs
      .map(function (tab, i) {
        return (
          '<div class="tab-pane fade' +
          (i === 0 ? " show active" : "") +
          '" id="tab-' +
          tab.id +
          '" role="tabpanel" aria-labelledby="tab-' +
          tab.id +
          '-btn" tabindex="0">' +
          tab.body +
          "</div>"
        );
      })
      .join("");
    return (
      '<ul class="nav nav-tabs flex-wrap result-tabs" role="tablist">' +
      nav +
      '</ul><div class="tab-content pt-3">' +
      panes +
      "</div>"
    );
  }

  function renderDashboard(root, payload) {
    const vcf = payload.vcf || {};
    const errors = payload.errors || [];
    root.innerHTML =
      '<div class="row g-3 mb-3">' +
      kpi(
        payload.ok ? "ok" : "failed",
        escapeHtml(payload.source_filename || ""),
        payload.ok ? "text-success" : "text-danger"
      ) +
      kpi(escapeHtml(vcf.sample_id || "—"), "sample") +
      kpi(Number(vcf.n_snps || 0).toLocaleString(), "SNPs parsed") +
      kpi(Number(vcf.n_non_snp_skipped || 0).toLocaleString(), "non-SNPs skipped") +
      "</div>" +
      '<div class="row g-3 mb-3">' +
      kpi(
        escapeHtml(vcf.lifted_to ? (vcf.assembly || "?") + " → " + vcf.lifted_to : vcf.assembly || "unknown"),
        vcf.lifted_to
          ? "assembly (lifted for AADR/hg19)"
          : "assembly",
        vcf.lifted_to ? "text-success" : vcf.assembly === "GRCh38" ? "text-warning" : ""
      ) +
      kpi(Number(vcf.n_lifted || 0).toLocaleString(), "sites lifted") +
      kpi(Number(vcf.n_unmapped || 0).toLocaleString(), "unmapped in liftover") +
      "</div>" +
      (vcf.assembly === "GRCh38" && !vcf.lifted_to
        ? '<div class="alert alert-warning">This VCF looks like GRCh38 (often GSA-24v3 / gtc2vcf). Haplogroups still match by rsID and hg38 positions. Place hg38ToHg19.over.chain.gz under data/references/liftover/ so autosomal AADR sites line up.</div>'
        : "") +
      (errors.length
        ? '<div class="alert alert-danger">' + errors.map(escapeHtml).join(" · ") + "</div>"
        : "") +
      resultTabs(payload);
    bindSnpList(root, payload);
    bindHgMarkerFilters(root);
  }

  async function loadSamples(select, emptyLabel) {
    const res = await fetch("/api/samples");
    const data = await res.json();
    select.innerHTML = "<option value=\"\">" + (emptyLabel || "Choose a bundled sample…") + "</option>";
    (data.samples || []).forEach(function (name) {
      const opt = document.createElement("option");
      opt.value = name;
      opt.textContent = name;
      select.appendChild(opt);
    });
  }

  function flagsFromForm(form) {
    return {
      hominin: form.querySelector('[name="hominin"]').checked,
      caste: form.querySelector('[name="caste"]').checked,
      populations: form.querySelector('[name="populations"]').checked,
      ancestry: form.querySelector('[name="ancestry"]').checked,
      haplogroups: !form.querySelector('[name="haplogroups"]') || form.querySelector('[name="haplogroups"]').checked,
      disease: !form.querySelector('[name="disease"]') || form.querySelector('[name="disease"]').checked,
    };
  }

  function wireBundledOrUpload(select, fileInput) {
    if (!select || !fileInput) return;
    select.addEventListener("change", function () {
      if (select.value) {
        fileInput.value = "";
      }
    });
    fileInput.addEventListener("change", function () {
      if (fileInput.files && fileInput.files.length) {
        select.value = "";
      }
    });
  }

  function formatElapsed(seconds) {
    const m = Math.floor(seconds / 60);
    const s = seconds % 60;
    return m + ":" + String(s).padStart(2, "0");
  }

  function buildAnalyzeStages(flags) {
    const stages = ["Reading genotype file", "Aligning SNPs to AADR panel"];
    const panelOn =
      flags.hominin || flags.caste || flags.populations || flags.ancestry;
    if (panelOn) {
      stages.push("Building or loading AADR frequency caches");
    }
    if (flags.populations) stages.push("Population mixture weights");
    if (flags.caste) stages.push("Community reference (caste)");
    if (flags.ancestry) stages.push("Ancestry models (3-source and 5-source)");
    if (flags.hominin) stages.push("Archaic hominin overlap");
    if (flags.haplogroups) stages.push("Y and mtDNA haplogroups");
    if (flags.disease) stages.push("Disease polygenic scores");
    stages.push("Additional trait markers");
    stages.push("Writing SNP catalog");
    return stages;
  }

  function createAnalyzeProgress(flags) {
    const panel = $("#progress-panel");
    const elapsedEl = $("#progress-elapsed");
    const bar = $("#progress-bar");
    const stagesEl = $("#progress-stages");
    if (!panel || !elapsedEl || !bar || !stagesEl) {
      return { stop: function () {} };
    }
    const stages = buildAnalyzeStages(flags || {});
    stagesEl.innerHTML = stages
      .map(function (label, i) {
        return (
          '<li class="' +
          (i === 0 ? "is-active" : "is-pending") +
          '" data-stage="' +
          i +
          '">' +
          escapeHtml(label) +
          "</li>"
        );
      })
      .join("");
    panel.classList.remove("d-none");
    bar.classList.add("progress-bar-animated", "progress-bar-striped");
    bar.style.width = "5%";
    bar.setAttribute("aria-valuenow", "5");

    const started = Date.now();
    let stageIndex = 0;

    function panelOnSlowStage(f) {
      return f.hominin || f.caste || f.populations || f.ancestry;
    }

    const stageIntervalMs = panelOnSlowStage(flags) ? 12000 : 4500;

    function renderStages() {
      const items = stagesEl.querySelectorAll("li");
      items.forEach(function (li, i) {
        li.className =
          i < stageIndex ? "is-done" : i === stageIndex ? "is-active" : "is-pending";
      });
      const pct = Math.min(
        95,
        Math.round(5 + ((stageIndex + 1) / stages.length) * 90)
      );
      bar.style.width = pct + "%";
      bar.setAttribute("aria-valuenow", String(pct));
    }

    const tick = window.setInterval(function () {
      const secs = Math.floor((Date.now() - started) / 1000);
      elapsedEl.textContent = "Analyzing… " + formatElapsed(secs);
      if (stageIndex < stages.length - 1) {
        const nextThreshold = (stageIndex + 1) * stageIntervalMs / 1000;
        if (secs >= nextThreshold) {
          stageIndex += 1;
          renderStages();
        }
      }
    }, 500);

    renderStages();

    return {
      stop: function (ok) {
        window.clearInterval(tick);
        stageIndex = stages.length - 1;
        renderStages();
        bar.classList.remove("progress-bar-animated", "progress-bar-striped");
        bar.style.width = "100%";
        bar.setAttribute("aria-valuenow", "100");
        elapsedEl.textContent = ok ? "Done." : "Stopped.";
        window.setTimeout(function () {
          panel.classList.add("d-none");
        }, ok ? 800 : 0);
      },
    };
  }

  function setupLive() {
    const form = $("#analyze-form");
    if (!form) return;
    const status = $("#status");
    const results = $("#results");
    const sampleSelect = $("#sample-select");
    const otherSelect = $("#other-select");
    const fileInput = $("#vcf-file");
    const otherInput = $("#other-file");
    const runBtn = $("#run-btn");
    wireBundledOrUpload(sampleSelect, fileInput);
    wireBundledOrUpload(otherSelect, otherInput);
    Promise.all([
      loadSamples(sampleSelect, "Choose a bundled sample…"),
      otherSelect ? loadSamples(otherSelect, "None — skip relatedness") : Promise.resolve(),
    ]).catch(function () {
      status.textContent = "Could not list bundled samples.";
    });

    form.addEventListener("submit", async function (event) {
      event.preventDefault();
      const sample = sampleSelect.value;
      const otherSample = otherSelect ? otherSelect.value : "";
      const flags = flagsFromForm(form);
      const primaryFile = fileInput.files && fileInput.files[0];
      const otherFile = otherInput && otherInput.files && otherInput.files[0];
      runBtn.disabled = true;
      status.className = "form-text text-secondary mb-0 mt-2";
      status.textContent = "";
      let progress = null;
      try {
        let res;
        if (primaryFile || otherFile) {
          if (!primaryFile && !sample) {
            status.className = "alert alert-danger py-2 px-3 mt-2 mb-0";
            status.textContent = "Choose a bundled sample or upload a SNP VCF / 23andMe raw .txt file.";
            return;
          }
          progress = createAnalyzeProgress(flags);
          const fd = new FormData();
          if (primaryFile) fd.append("vcf", primaryFile);
          if (sample) fd.append("sample", sample);
          if (otherFile) fd.append("other", otherFile);
          if (otherSample) fd.append("other_sample", otherSample);
          fd.append("hominin", String(flags.hominin));
          fd.append("caste", String(flags.caste));
          fd.append("populations", String(flags.populations));
          fd.append("ancestry", String(flags.ancestry));
          fd.append("haplogroups", String(flags.haplogroups));
          fd.append("disease", String(flags.disease));
          res = await fetch("/api/analyze", { method: "POST", body: fd });
        } else if (sample) {
          progress = createAnalyzeProgress(flags);
          res = await fetch("/api/analyze-sample", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(Object.assign({
              filename: sample,
              other_filename: otherSample || "",
            }, flags)),
          });
        } else {
          status.className = "alert alert-danger py-2 px-3 mt-2 mb-0";
          status.textContent = "Choose a bundled sample or upload a SNP VCF / 23andMe raw .txt file.";
          return;
        }
        const payload = await res.json();
        if (!res.ok && !payload.source_filename) {
          throw new Error(payload.error || "Analyze failed");
        }
        if (progress) progress.stop(true);
        renderDashboard(results, payload);
        results.classList.remove("d-none");
        status.textContent = payload.ok ? "Analysis complete." : "Finished with errors.";
        status.className = payload.ok
          ? "form-text text-success mb-0 mt-2"
          : "alert alert-danger py-2 px-3 mt-2 mb-0";
        if (payload.timings && payload.timings.total_ms) {
          status.textContent +=
            " Total " + Math.round(payload.timings.total_ms / 1000) + "s.";
        }
      } catch (err) {
        if (progress) progress.stop(false);
        status.className = "alert alert-danger py-2 px-3 mt-2 mb-0";
        status.textContent = err.message || String(err);
      } finally {
        runBtn.disabled = false;
      }
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    const results = $("#results");
    if (window.ANALYSIS_PAYLOAD && results) {
      $("#live-controls") && $("#live-controls").classList.add("d-none");
      renderDashboard(results, window.ANALYSIS_PAYLOAD);
      results.classList.remove("d-none");
      return;
    }
    setupLive();
  });

  window.renderDnaDashboard = renderDashboard;
})();
