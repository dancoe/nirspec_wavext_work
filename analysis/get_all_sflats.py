import crds
import crds.core.log
import crds.jwst.locate

ctx = crds.get_default_context()
rmap = crds.get_cached_mapping(ctx)
nirspec_imap = rmap.get_imap("nirspec")
sflat_rmap = nirspec_imap.get_rmap("sflat")

sflats = {}

# We need to find all combinations of EXP_TYPE, DETECTOR, GRATING, FILTER
# that map to an sflat in this rmap.
# Using sflat_rmap.reference_names() gives all reference filenames in the rmap.
all_refs = sflat_rmap.reference_names()

print("Found", len(all_refs), "sflat references in", ctx)
with open("/tmp/crds_sflats.txt", "w") as f:
    for ref in all_refs:
        # get parkeys from the filename by querying crds
        try:
            headers = crds.get_file_properties("jwst", ref)
            print(f"{ref}: {headers}", file=f)
        except Exception as e:
            # Maybe just try reading header if it's cached, otherwise download?
            pass

