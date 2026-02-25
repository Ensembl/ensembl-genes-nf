# Translon Consensus Pipeline Assets

This directory contains templates and static assets for the pipeline.

## Files

### report_template.html
Jinja2 HTML template for generating the final consensus report.

**Features:**
- Interactive tables with sorting/filtering (DataTables.js)
- Interactive charts (Plotly.js)
- Responsive design
- Self-contained (uses CDN for JS/CSS)

**Template Variables:**
- `run_name`: Pipeline run identifier
- `timestamp`: Report generation timestamp
- `ucsc_session_url`: UCSC Genome Browser session URL
- `total_samples`: Total number of samples
- `multi_tool_samples`: Number of samples with 2+ tools
- `single_tool_samples`: Number of samples with only 1 tool
- `unique_tools`: Number of unique tools across all samples
- `all_tools`: List of all tool names
- `tool_names`: Tool names for charts
- `tool_counts`: Tool usage counts
- `samples_with_consensus`: List of samples with consensus analysis
- `skipped_samples`: List of single-tool samples (no consensus)

**Customization:**
To customize the report appearance, edit this template file. The template uses:
- Jinja2 syntax for templating
- CSS Grid for layout
- DataTables.js for interactive tables
- Plotly.js for charts

**Usage:**
The template is automatically loaded by `generate_html_report.py`. You can also specify a custom template with the `-t` flag:

```bash
generate_html_report.py -i results/ -o report.html -t custom_template.html
```
