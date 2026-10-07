import streamlit as st
import subprocess
import tempfile
import pathlib
import zipfile
import itertools
import gzip
import io


st.set_page_config(page_title="SequenceDot", layout="centered")
st.title("SequenceDot")
st.write("K-mer based dotplot generator for biological sequences.")

# ── Mode ──────────────────────────────────────────────────────────────────
mode = st.radio(
    "Comparison mode",
    ["Two sequences", "Batch (multi-FASTA)"],
    horizontal=True
)

# ── File upload ───────────────────────────────────────────────────────────
if mode == "Two sequences":
    seq1_file = st.file_uploader(
        "Sequence 1 (.fasta, .fa, .fna, or .gz)",
        type=["fasta", "fa", "fna", "gz"]
    )
    seq2_file = st.file_uploader(
        "Sequence 2 (.fasta, .fa, .fna, or .gz)",
        type=["fasta", "fa", "fna", "gz"]
    )
    multi_file     = None
    include_self   = False
else:
    multi_file = st.file_uploader(
        "Multi-sequence FASTA (.fasta, .fa, .fna, or .gz). For large datasets, use the command-line tool.",
        type=["fasta", "fa", "fna", "gz"]
    )
    include_self = st.checkbox("Include self-comparisons", value=False)
    seq1_file = seq2_file = None

# ── Parameters ────────────────────────────────────────────────────────────
kmer       = st.slider("K-mer size", min_value=1,   max_value=100,  value=11)
point_size = st.slider("Dot size",   min_value=0.1, max_value=20.0, value=1.0, step=0.1)
alphabet   = st.radio("Alphabet",        ["DNA", "RNA", "AA"],          horizontal=True)
strand     = st.radio("Strand",          ["forward", "reverse", "both"], horizontal=True)
fmt        = st.radio("Download format", ["PNG", "PDF", "SVG"],          horizontal=True)

# ── Session state ─────────────────────────────────────────────────────────
for key, default in [("plots", []), ("plot_index", 0)]:
    if key not in st.session_state:
        st.session_state[key] = default

# ── Helpers ───────────────────────────────────────────────────────────────

def parse_fasta(file_obj):
    """Return list of (seq_id, sequence) from an UploadedFile."""
    raw = file_obj.getvalue()
    if file_obj.name.endswith(".gz"):
        raw = gzip.decompress(raw)
    content = raw.decode("utf-8", errors="replace")

    records, current_id, current_seq = [], None, []
    for line in content.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, "".join(current_seq)))
            current_id = line[1:].split()[0]
            current_seq = []
        elif line:
            current_seq.append(line)
    if current_id is not None:
        records.append((current_id, "".join(current_seq)))
    return records


def safe_name(name):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)


def write_fasta(directory, seq_id, sequence):
    path = pathlib.Path(directory) / f"{safe_name(seq_id)}.fasta"
    path.write_text(f">{seq_id}\n{sequence}\n")
    return path


def run_pair(seq1_path, seq2_path, out_path, kmer, alphabet, strand, point_size):
    return subprocess.run(
        [
            "seqdot",
            str(seq1_path), str(seq2_path),
            "--kmer",       str(kmer),
            "--alphabet",   alphabet,
            "--strand",     strand,
            "--point-size", str(point_size),
            "--output",     str(out_path),
        ],
        capture_output=True, text=True
    )


def make_pair_plot(pair_dir, id1, seq1, id2, seq2, ext, kmer, alphabet, strand, point_size):
    """Generate one plot; return (output_path, error_string)."""
    p1 = write_fasta(pair_dir, id1, seq1)
    p2 = write_fasta(pair_dir, id2, seq2)
    out = pathlib.Path(pair_dir) / f"{safe_name(id1)}_vs_{safe_name(id2)}.{ext}"
    result = run_pair(p1, p2, out, kmer, alphabet, strand, point_size)
    if result.returncode != 0:
        return None, result.stderr or result.stdout
    return (out if out.exists() else None), ""


# ── Generate ──────────────────────────────────────────────────────────────

if st.button("Generate dotplot(s)", type="primary"):
    st.session_state.plots      = []
    st.session_state.plot_index = 0

    fmt_ext  = fmt.lower()
    mime_map = {"png": "image/png", "pdf": "application/pdf", "svg": "image/svg+xml"}
    plots    = []
    errors   = []

    if mode == "Two sequences":
        if seq1_file is None or seq2_file is None:
            st.error("Please upload both sequence files.")
        else:
            d = tempfile.mkdtemp()
            p1 = pathlib.Path(d) / seq1_file.name
            p2 = pathlib.Path(d) / seq2_file.name
            p1.write_bytes(seq1_file.getvalue())
            p2.write_bytes(seq2_file.getvalue())

            dl_out = pathlib.Path(d) / f"dotplot.{fmt_ext}"
            result = run_pair(p1, p2, dl_out, kmer, alphabet, strand, point_size)

            if result.returncode != 0:
                st.error(f"SequenceDot error:\n{result.stderr or result.stdout}")
            else:
                if fmt_ext == "png":
                    png_out = dl_out
                else:
                    png_out = pathlib.Path(d) / "dotplot_display.png"
                    run_pair(p1, p2, png_out, kmer, alphabet, strand, point_size)

                if dl_out.exists() and png_out.exists():
                    plots.append({
                        "name":    f"{seq1_file.name} vs {seq2_file.name}",
                        "png":     str(png_out),
                        "dl":      str(dl_out),
                        "dl_ext":  fmt_ext,
                        "dl_mime": mime_map[fmt_ext],
                    })

    else:  # Batch
        if multi_file is None:
            st.error("Please upload a multi-sequence FASTA file.")
        else:
            records = parse_fasta(multi_file)

            if len(records) < 2:
                st.error("The file must contain at least 2 sequences.")
            else:
                if include_self:
                    pairs = [
                        (records[i], records[j])
                        for i in range(len(records))
                        for j in range(i, len(records))
                    ]
                else:
                    pairs = list(itertools.combinations(records, 2))


                bar = st.progress(0, text="Generating plots…")

                for i, ((id1, seq1), (id2, seq2)) in enumerate(pairs):
                    bar.progress(
                        (i + 1) / len(pairs),
                        text=f"Comparing {id1} vs {id2}…"
                    )
                    d = tempfile.mkdtemp()

                    dl_path, err = make_pair_plot(
                        d, id1, seq1, id2, seq2,
                        fmt_ext, kmer, alphabet, strand, point_size
                    )
                    if dl_path is None:
                        errors.append(f"{id1} vs {id2}: {err}")
                        continue

                    if fmt_ext == "png":
                        png_path = dl_path
                    else:
                        png_path, _ = make_pair_plot(
                            d, id1, seq1, id2, seq2,
                            "png", kmer, alphabet, strand, point_size
                        )

                    if png_path and png_path.exists():
                        plots.append({
                            "name":    f"{id1} vs {id2}",
                            "png":     str(png_path),
                            "dl":      str(dl_path),
                            "dl_ext":  fmt_ext,
                            "dl_mime": mime_map[fmt_ext],
                        })

                bar.empty()
                for e in errors:
                    st.warning(e)

    st.session_state.plots = plots

# ── Display ───────────────────────────────────────────────────────────────

if st.session_state.plots:
    plots = st.session_state.plots
    n     = len(plots)
    idx   = st.session_state.plot_index
    curr  = plots[idx]

    # Navigation row
    c_prev, c_label, c_next = st.columns([1, 6, 1])
    with c_prev:
        if st.button("◀", disabled=(idx == 0), use_container_width=True):
            st.session_state.plot_index -= 1
            st.rerun()
    with c_label:
        st.markdown(
            f"<p style='text-align:center;margin-top:6px'>"
            f"<b>{curr['name']}</b> &nbsp;·&nbsp; {idx + 1} / {n}</p>",
            unsafe_allow_html=True
        )
    with c_next:
        if st.button("▶", disabled=(idx == n - 1), use_container_width=True):
            st.session_state.plot_index += 1
            st.rerun()

    # Plot
    st.image(curr["png"], use_container_width=True)

    # Individual download
    dl_path = pathlib.Path(curr["dl"])
    if dl_path.exists():
        st.download_button(
            label=f"⬇ Download this plot ({curr['dl_ext'].upper()})",
            data=dl_path.read_bytes(),
            file_name=dl_path.name,
            mime=curr["dl_mime"],
            use_container_width=True
        )

    # Download all as ZIP
    if n > 1:
        st.divider()
        zip_buf = io.BytesIO()
        with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for p in plots:
                dp = pathlib.Path(p["dl"])
                if dp.exists():
                    zf.write(dp, arcname=dp.name)
        zip_buf.seek(0)
        st.download_button(
            label=f"⬇ Download all {n} plots as ZIP ({curr['dl_ext'].upper()})",
            data=zip_buf.getvalue(),
            file_name="sequencedot_plots.zip",
            mime="application/zip",
            use_container_width=True
        )
