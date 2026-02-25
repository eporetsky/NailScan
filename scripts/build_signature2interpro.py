#!/usr/bin/env python3
"""
Build signature -> InterPro ID mapping from interpro.xml.gz.
Reads from stdin or first arg path; writes db\tsignature\tinterpro_id to stdout.
Used by download_interpro_mappings.sh to create data/signature2interpro.tsv.
"""
import re
import sys
import gzip

# Our config keys -> InterPro XML member db names (as in <db_xref db="...">)
DB_TO_XML = {
    "pfam": "PFAM",
    "ncbifam": "NCBIFAM",
    "superfamily": "SSF",
    "hamap": "HAMAP",
    "sfld": "SFLD",
    "pirsf": "PIRSF",
    "pirsr": "PIRSR",
    "panther": "PANTHER",
    "cath": "CATHGENE3D",
    "antifam": "ANTIFAM",
}
XML_DBS = set(DB_TO_XML.values())


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if path and path != "-":
        f = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        f = sys.stdin

    interpro_re = re.compile(r'<interpro\s+id="(IPR[^"]+)"')
    member_re1 = re.compile(r'<db_xref\s+[^>]*\bdb="([^"]+)"[^>]*\bdbkey="([^"]+)"')
    member_re2 = re.compile(r'<db_xref\s+[^>]*\bdbkey="([^"]+)"[^>]*\bdb="([^"]+)"')
    current_ipr = None

    for line in f:
        m = interpro_re.search(line)
        if m:
            current_ipr = m.group(1)
            continue
        if current_ipr and "<member_list>" in line:
            continue
        if current_ipr and "</member_list>" in line:
            current_ipr = None
            continue
        if current_ipr:
            mm = member_re1.search(line)
            if mm:
                db, dbkey = mm.group(1), mm.group(2)
            else:
                mm = member_re2.search(line)
                if mm:
                    dbkey, db = mm.group(1), mm.group(2)
                else:
                    mm = None
            if mm and db in XML_DBS and dbkey:
                print(db, dbkey, current_ipr, sep="\t")

    if path and path != "-":
        f.close()


if __name__ == "__main__":
    main()
