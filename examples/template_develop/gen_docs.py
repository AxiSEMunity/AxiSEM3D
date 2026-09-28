#!/usr/bin/env python3
#
#  gen_docs.py
#  AxiSEM3D — generate human-readable parameter docs from annotated inparam YAML.
#
#  Safe to run from any directory. Scope: only input/inparam*.yaml next to this
#  script (examples/template_develop/input/). Writes MyST-friendly markdown to
#  doc/sphinx/source/parameter/<basename>.md (see OUTPUT_DIR).

import glob
import os
import re

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_GLOB = os.path.join(SCRIPT_DIR, "input", "inparam*.yaml")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "../../doc/sphinx/source/parameter")

# Section banner: ##### title #####
SECTION_LINE = re.compile(r"^#{7,}\s*(.+?)\s*#{7,}\s*$")
# YAML key: optional value (leaf if value non-empty after strip)
PARAM_LINE = re.compile(r"^(\s*)([A-Za-z_][\w.]*)\s*:\s*(.*)$")
META_START = re.compile(r"(?i)^(what|type|only|note)\s*:\s*(.*)$")
SECTION_BLURB = re.compile(r"(?i)^#\s*parameters\s+for\s+(.+)$")


def collect_comment_block_above(lines, param_idx):
    """Lines are '# ...' metadata immediately above lines[param_idx]."""
    block = []
    j = param_idx - 1
    while j >= 0:
        raw = lines[j]
        stripped = raw.strip()
        if not stripped:
            j -= 1
            continue
        if stripped.startswith("#"):
            block.append(raw)
            j -= 1
            continue
        m = PARAM_LINE.match(raw)
        if m:
            _ind, _key, val = m.groups()
            if val.strip() == "":
                break
            break
        break
    block.reverse()
    return block


def parse_meta(comment_lines):
    what = typ = only = ""
    note_chunks = []
    in_note = False

    for raw in comment_lines:
        stripped = raw.strip()
        if not stripped.startswith("#"):
            continue
        inner = stripped[1:].strip()
        mm = META_START.match(inner)
        if mm:
            field, rest = mm.group(1).lower(), mm.group(2).strip()
            in_note = field == "note"
            if field == "what":
                what = rest
            elif field == "type":
                typ = rest
            elif field == "only":
                only = rest
            elif field == "note":
                note_chunks.append(rest)
            continue
        if in_note or note_chunks:
            note_chunks.append(inner)
            in_note = True
            continue
    note = "\n".join(note_chunks).strip() if note_chunks else ""
    return what, typ, only, note


def sniff_subsection_preamble(lines, block_idx):
    """Free-form # lines directly above a `key:` block header (empty value)."""
    j = block_idx - 1
    parts = []
    while j >= 0:
        stripped = lines[j].strip()
        if not stripped:
            j -= 1
            continue
        if stripped.startswith("#"):
            inner = stripped[1:].strip()
            if META_START.match(inner):
                break
            if re.match(r"(?i)^parameters\s+for\s+", inner):
                break
            parts.append(inner)
            j -= 1
            continue
        break
    parts.reverse()
    return "\n".join(parts).strip() if parts else ""


def section_blurb_from_line(line):
    m = SECTION_BLURB.match(line.strip())
    if m:
        return m.group(1).strip()
    return None


def extract_section_title(banner_line):
    m = SECTION_LINE.match(banner_line.strip())
    return m.group(1).strip() if m else None


def yaml_leaf_value(value_part):
    return value_part.strip() != ""


def md_default_value(val):
    """Inline or fenced literal for the YAML RHS (no # default: in source files)."""
    if not val.strip():
        return ""
    if "\n" in val or "`" in val:
        return "\n\n```\n" + val.rstrip() + "\n```\n"
    return f" `{val}`"


def format_parameter_block(key, what, typ, only, note, yaml_value, in_subsection):
    """Use ### under ##; use #### only after a ### subsection (consecutive MyST levels)."""
    level = "####" if in_subsection else "###"
    parts_out = [f"{level} `{key}`\n\n"]
    if what:
        parts_out.append(f"**What:** {what}\n\n")
    if typ:
        parts_out.append(f"**Type:** {typ}\n\n")
    if only:
        parts_out.append(f"**Only:** {only}\n\n")
    parts_out.append("**Default:**")
    parts_out.append(md_default_value(yaml_value))
    parts_out.append("\n\n")
    if note:
        parts_out.append("**Note:**\n\n")
        parts_out.append(note + "\n\n")
    parts_out.append("\n")
    return "".join(parts_out)


def document_inparam_yaml(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    section_title = None
    section_blurb = None
    # stack of (indent, key, preamble) for YAML block headers with empty value
    parent_stack = []
    out = []
    last_section_emitted = None
    last_subsection_signature = None

    for i, line in enumerate(lines):
        banner = extract_section_title(line)
        if banner:
            section_title = banner
            section_blurb = None
            parent_stack.clear()
            continue

        blurb = section_blurb_from_line(line)
        if blurb is not None:
            section_blurb = blurb
            continue

        m = PARAM_LINE.match(line)
        if not m:
            continue

        indent_s, key, val_part = m.groups()
        leaf_indent = len(indent_s.expandtabs())

        if not yaml_leaf_value(val_part):
            while parent_stack and parent_stack[-1][0] >= leaf_indent:
                parent_stack.pop()
            pre = sniff_subsection_preamble(lines, i)
            parent_stack.append((leaf_indent, key, pre))
            continue

        while parent_stack and parent_stack[-1][0] >= leaf_indent:
            parent_stack.pop()

        comments = collect_comment_block_above(lines, i)
        what, typ, only, note = parse_meta(comments)
        yaml_value = val_part.strip()

        if section_title and section_title != last_section_emitted:
            out.append(f"## {section_title}\n\n")
            if section_blurb:
                out.append(f"Parameters for {section_blurb}.\n\n")
            last_section_emitted = section_title
            last_subsection_signature = None

        sig = None
        if len(parent_stack) >= 2:
            sig = ("depth", parent_stack[-1][0], parent_stack[-1][1])
        elif parent_stack and parent_stack[-1][2]:
            sig = ("preamble", parent_stack[-1][0], parent_stack[-1][1])

        if sig and sig != last_subsection_signature:
            sub_key = parent_stack[-1][1]
            sub_pre = parent_stack[-1][2]
            out.append(f"### {sub_key}\n\n")
            if sub_pre:
                out.append(f"{sub_pre}\n\n")
            last_subsection_signature = sig

        in_subsection = sig is not None
        out.append(
            format_parameter_block(
                key, what, typ, only, note, yaml_value, in_subsection
            )
        )

    return "".join(out)

def writeExampleCustomName(path,name):
        base = os.path.basename(name)
        base_cleaned = base.replace(" ", "_")
        stem = os.path.splitext(base_cleaned)[0]
        out_path = os.path.join(OUTPUT_DIR, stem + ".md")
        body = document_inparam_yaml(path)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(f"# {base}\n\n")
            f.write(body)
        print(f"Wrote {out_path}")

def main():
    os.chdir(SCRIPT_DIR)
    paths = sorted(glob.glob(INPUT_GLOB))
    if not paths:
        raise SystemExit("No inparam*.yaml files under input/")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for path in paths:
        base = os.path.basename(path)
        stem = os.path.splitext(base)[0]
        stem = stem.replace("inparam.", "")
        out_path = os.path.join(OUTPUT_DIR, stem + ".md")
        body = document_inparam_yaml(path)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(f"# {base}\n\n")
            f.write(body)
        print(f"Wrote {out_path}")

    # Stations, sources, models examples should not necessarily 
    # be named after their example files.
    writeExampleCustomName(os.path.join(SCRIPT_DIR, "input", "ex06_list_of_station_groups.yaml"),"List of Stations")
    writeExampleCustomName(os.path.join(SCRIPT_DIR, "input", "ex05c_list_of_sources.yaml"),"List of Sources")
    writeExampleCustomName(os.path.join(SCRIPT_DIR, "input", "ex06_list_of_element_groups.yaml"),"List of Element Groups")
    writeExampleCustomName(os.path.join(SCRIPT_DIR, "input", "ex05c_list_of_3D_models.yaml"),"List of 3D Models")

if __name__ == "__main__":
    main()
