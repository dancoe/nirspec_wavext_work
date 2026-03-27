from astroquery.mast import Observations
import pandas as pd
import numpy as np

prog_id = '1492'
print(f"Querying PID {prog_id} NIRSPEC/SLIT observations...")

obs_table = Observations.query_criteria(
    proposal_id=prog_id,
    instrument_name='NIRSPEC/SLIT',
    intentType='science'
)

df = obs_table.to_pandas()
print(f"Total rows found: {len(df)}")

# Summary for Level 3
l3_df = df[df['calib_level'] == 3].copy()
l3_df['grating'] = l3_df['filters'].str.split(';').str[-1]
l3_summary = l3_df.groupby('grating').agg(
    count=('obs_id', 'count'),
    t_exptime=('t_exptime', 'mean')
).reset_index()

# Summary for Level 2
l2_df = df[df['calib_level'] == 2].copy()
l2_df['grating'] = l2_df['filters'].str.split(';').str[-1]
l2_summary = l2_df.groupby('grating').agg(
    count=('obs_id', 'count'),
    t_exptime=('t_exptime', 'mean')
).reset_index()

output_file = '/Users/dcoe/NIRSpec/wavext/pipeline/notes/DATA.md'

with open(output_file, 'w') as f:
    f.write(f"# DATA.md - Program ID {prog_id} NIRSpec FS Observations\n\n")
    f.write(f"**Target:** {df.iloc[0]['target_name']}\n\n")
    
    f.write("## Usage Strategy\n\n")
    f.write("1.  **Calibration Target Selection:** Focus on the high-SNR sources for each grating.\n")
    f.write("2.  **Extended Extraction:** Process these observations through the modified pipeline with the `wavelengthrange_extended.asdf` override.\n")
    f.write("3.  **Cross-Check:** Compare the extended cutouts on NRS2 with the nominal NRS1 data for G140/G235/G395 combinations to measure the k(λ), a(λ), and b(λ) functions.\n\n")

    f.write("## High-Level (Level 3) Products Summary\n\n")
    f.write("| Grating | Count | Avg Exp Time (s) |\n")
    f.write("| --- | --- | --- |\n")
    for i, row in l3_summary.iterrows():
        f.write(f"| {row['grating']} | {row['count']} | {row['t_exptime']:.2f} |\n")

    f.write("\n## Level 2 Files Summary\n\n")
    f.write("| Grating | Count | Avg Exp Time (s) |\n")
    f.write("| --- | --- | --- |\n")
    for i, row in l2_summary.iterrows():
        f.write(f"| {row['grating']} | {row['count']} | {row['t_exptime']:.2f} |\n")

    f.write("\n\n## Complete Observations List\n\n")
    # Write as a table
    table_cols = ['obs_id', 'filters', 't_exptime', 'calib_level']
    header = "| " + " | ".join([c.replace('_', ' ').title() for c in table_cols]) + " |\n"
    divider = "| " + " | ".join(["---"] * len(table_cols)) + " |\n"
    f.write(header)
    f.write(divider)
    
    # Sort for display
    df_disp = df.sort_values(by=['calib_level', 'filters', 't_exptime'], ascending=[False, True, True])
    for i, row in df_disp.iterrows():
        vals = [str(row[c]) for c in table_cols]
        f.write("| " + " | ".join(vals) + " |\n")

print(f"Notes saved to {output_file}")
