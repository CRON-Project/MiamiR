
build_app_help <- function(
    pkg_root = ".",
    out_dir  = file.path(pkg_root, "inst", "app_help"),
    pkg_name = "MiamiR",
    verbose  = TRUE
) {
  man_dir <- file.path(pkg_root, "man")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  if (!length(rd_files)) {
    message("No .Rd files found in: ", man_dir)
    return(invisible(out_dir))
  }

  topics <- sub("\\.Rd$", "", basename(rd_files))

  # 1) Write a nice stylesheet the app WILL have
  css_path <- file.path(out_dir, "style.css")
  writeLines(c(
    "/* MiamiR app help styling */",
    ":root { --bg:#0b0f19; --card:#111827; --text:#e5e7eb; --muted:#9ca3af; --link:#60a5fa; --border:#1f2937; }",
    "html, body { height: 100%; }",
    "body { margin:0; font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;",
    "       background: var(--bg); color: var(--text); line-height: 1.55; }",
    ".topbar { position: sticky; top:0; z-index:10; background: rgba(17,24,39,0.95);",
    "          border-bottom: 1px solid var(--border); padding: 12px 16px; display:flex; gap:12px; align-items:center; }",
    ".topbar a { color: var(--text); text-decoration:none; font-weight: 650; }",
    ".topbar .crumb { color: var(--muted); font-weight: 500; }",
    ".container { max-width: 980px; margin: 18px auto; padding: 0 16px 28px; }",
    ".card { background: var(--card); border: 1px solid var(--border); border-radius: 14px; padding: 18px; box-shadow: 0 10px 30px rgba(0,0,0,0.25); }",
    "a { color: var(--link); }",
    "h1,h2,h3 { line-height: 1.15; margin-top: 1.1em; }",
    "pre, code { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, \"Liberation Mono\", \"Courier New\", monospace; }",
    "pre { background: #0a1220; border: 1px solid var(--border); padding: 12px; border-radius: 12px; overflow:auto; }",
    "code { background: rgba(96,165,250,0.12); padding: 0.1em 0.35em; border-radius: 8px; }",
    "table { border-collapse: collapse; width: 100%; }",
    "th, td { border: 1px solid var(--border); padding: 8px 10px; vertical-align: top; }",
    ".search { display:flex; gap:10px; margin: 14px 0 0; }",
    ".search input { width: 100%; padding: 10px 12px; border-radius: 12px; border: 1px solid var(--border);",
    "               background: #0a1220; color: var(--text); outline: none; }",
    ".list a { display:block; padding: 8px 10px; border-radius: 10px; text-decoration:none; }",
    ".list a:hover { background: rgba(96,165,250,0.12); }",
    ".muted { color: var(--muted); }"
  ), css_path, useBytes = TRUE)

  postprocess <- function(html_file, topic) {
    x <- readLines(html_file, warn = FALSE, encoding = "UTF-8")

    # Swap the default Rd2HTML stylesheet reference (usually R.css) to our shipped CSS
    x <- gsub('href="R.css"', 'href="style.css"', x, fixed = TRUE)

    # Ensure we *have* our stylesheet link even if Rd2HTML output differs
    if (!any(grepl('href="style.css"', x, fixed = TRUE))) {
      head_i <- which(grepl("<head>", x, fixed = TRUE))
      if (length(head_i)) {
        x <- append(x, values = '  <link rel="stylesheet" type="text/css" href="style.css">', after = head_i[1])
      }
    }

    # Wrap content in a nice layout + header
    body_i <- which(grepl("<body>", x, fixed = TRUE))
    if (length(body_i)) {
      header <- c(
        "<body>",
        sprintf('<div class="topbar"><a href="index.html">%s Help</a><span class="crumb">/ %s</span></div>', pkg_name, topic),
        '<div class="container"><div class="card">'
      )
      # Replace the first <body> line with our header block
      x[body_i[1]] <- header[1]
      x <- append(x, header[-1], after = body_i[1])

      end_i <- which(grepl("</body>", x, fixed = TRUE))
      if (length(end_i)) {
        # insert closing divs BEFORE </body>
        x <- append(x, c("</div></div>"), after = end_i[1] - 1L)
      } else {
        x <- c(x, "</div></div>")
      }
    }

    writeLines(x, html_file, useBytes = TRUE)
  }


  # 2) Convert + progress messages
  n <- length(rd_files)
  if (verbose) message("Building HTML help (", n, " topics) -> ", out_dir)

  for (i in seq_along(rd_files)) {
    topic <- topics[i]
    out_file <- file.path(out_dir, paste0(topic, ".html"))

    if (verbose) message(sprintf("[%d/%d] Converting: %s", i, n, topic))

    rd <- tools::parse_Rd(rd_files[i])
    tools::Rd2HTML(rd, out = out_file)

    postprocess(out_file, topic)
  }

  # 3) Create a nice index.html with search
  idx <- file.path(out_dir, "index.html")
  links <- paste0('<a href="', topics, '.html">', topics, "</a>", collapse = "\n")

  writeLines(c(
    "<!doctype html>",
    "<html><head>",
    '  <meta charset="utf-8">',
    sprintf("  <title>%s Help</title>", pkg_name),
    '  <link rel="stylesheet" type="text/css" href="style.css">',
    "</head><body>",
    sprintf('<div class="topbar"><a href="index.html">%s Help</a><span class="crumb muted">Reference</span></div>', pkg_name),
    '<div class="container">',
    '<div class="card">',
    sprintf("<h1>%s documentation</h1>", pkg_name),
    '<p class="muted">Click a topic or search below.</p>',
    '<div class="search"><input id="q" type="text" placeholder="Search topics..."></div>',
    '<div id="list" class="list">',
    links,
    "</div>",
    "</div></div>",
    "<script>",
    "  const q = document.getElementById('q');",
    "  const list = document.getElementById('list');",
    "  const items = Array.from(list.querySelectorAll('a'));",
    "  q.addEventListener('input', () => {",
    "    const s = q.value.toLowerCase();",
    "    items.forEach(a => {",
    "      a.style.display = a.textContent.toLowerCase().includes(s) ? '' : 'none';",
    "    });",
    "  });",
    "</script>",
    "</body></html>"
  ), idx, useBytes = TRUE)

  message("Wrote HTML help to: ", out_dir)
  invisible(out_dir)
}

# run it
build_app_help()

