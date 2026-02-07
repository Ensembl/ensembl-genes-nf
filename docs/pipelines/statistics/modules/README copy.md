# Statistics Pipeline Module Documentation

Welcome to the comprehensive module documentation for the Ensembl Statistics Pipeline!

## 🚀 Quick Start

### First Time Here?

1. **Start with the Index**: [index.md](index.md) provides an overview of all 14 modules
2. **Pick a module**: Choose based on your needs (see categories below)
3. **Read the documentation**: Each module has detailed docs with examples
4. **Try the examples**: Copy-paste code and test locally

### Looking for Something Specific?

| I want to... | Go to... |
|--------------|----------|
| Understand what each module does | [Module Index - Overview Section](index.md#-module-overview) |
| See how modules work together | [Module Index - Module Usage by Workflow](index.md#-module-usage-by-workflow) |
| Learn about a specific module | [Module Index - Alphabetical List](index.md#alphabetical) |
| Troubleshoot an error | Specific module doc → Error Handling section |
| Optimize performance | Specific module doc → Performance Considerations |
| Write tests | Specific module doc → Testing section |

## 📚 Available Documentation

### Complete Module Documentation

| Module | Purpose | Documentation | Status |
|--------|---------|---------------|--------|
| **DB_METADATA** | Extract metadata from Ensembl databases | [View Docs](db-metadata.md) | ✅ Complete (518 lines) |
| **FETCH_GENOME** | Download genome assemblies from NCBI | [View Docs](fetch-genome.md) | ✅ Complete (642 lines) |
| **CLEANING** | Clean up temporary files | [View Docs](cleaning.md) | ✅ Complete (495 lines) |

### Index Documentation

| Document | Purpose | Lines |
|----------|---------|-------|
| **index.md** | Overview of all modules, navigation hub | 243 lines |

**Total Documentation**: 1,898 lines across 4 files

## 🗂️ Module Categories

### Data Retrieval Modules

Modules that fetch input data from databases and external sources:

- **[DB_METADATA](db-metadata.md)** ✅ - Extract metadata (taxon_id, GCA, production_name)
- **FETCH_GENOME** ✅ - Download genomes from NCBI
- **FETCH_PROTEINS** 🔲 - Extract protein sequences

### Analysis Modules

Modules that perform quality assessment and statistics:

- **BUSCO_DATASET** 🔲 - Select BUSCO lineage
- **BUSCO_GENOME_LINEAGE** 🔲 - Run BUSCO on genomes
- **BUSCO_PROTEIN_LINEAGE** 🔲 - Run BUSCO on proteins
- **OMAMER_HOG** 🔲 - Generate HOG assignments
- **OMARK** 🔲 - Assess proteome quality
- **RUN_STATISTICS** 🔲 - Generate assembly statistics
- **RUN_ENSEMBL_META** 🔲 - Generate metadata statistics

### Database Integration Modules

Modules for storing results in databases:

- **BUSCO_CORE_METAKEYS** 🔲 - Insert BUSCO results
- **POPULATE_DB** 🔲 - Load statistics

### Utility Modules

Supporting modules:

- **[CLEANING](cleaning.md)** ✅ - Cleanup utilities

**Legend**: ✅ Complete documentation | 🔲 Coming soon

## 📖 Documentation Structure

Each complete module documentation includes:

### 1. Overview
- Purpose and functionality
- Module location
- Quick summary

### 2. Core Documentation
- **Inputs**: Channel inputs, metadata maps, parameters
- **Outputs**: Channel outputs, file formats, locations
- **Process Configuration**: Resources, containers, directives

### 3. Implementation
- Script logic and tool usage
- Caching strategies
- Data handling patterns

### 4. Practical Guides
- **Usage Examples**: Complete workflow examples
- **Error Handling**: Common errors with solutions
- **Best Practices**: Recommended patterns
- **Performance**: Execution times, optimization

### 5. Testing & Troubleshooting
- Unit tests with test data
- Debug procedures
- Validation steps

### 6. Integration
- Upstream/downstream dependencies
- Data flow diagrams
- Related modules

## 🎯 Common Use Cases

### Use Case 1: Understanding Module Inputs

**Scenario**: "What inputs does FETCH_GENOME need?"

**Solution**:
1. Open [fetch-genome.md](fetch-genome.md)
2. Go to **Inputs** section
3. See complete metadata map with all fields
4. Check **Parameters** table for configuration options

### Use Case 2: Troubleshooting Errors

**Scenario**: "FETCH_GENOME failed with 'Connection timeout'"

**Solution**:
1. Open [fetch-genome.md](fetch-genome.md)
2. Go to **Error Handling** → **Network Timeout**
3. Follow suggested solutions
4. If needed, check **Troubleshooting** section

### Use Case 3: Optimizing Performance

**Scenario**: "Can I run multiple FETCH_GENOME processes in parallel?"

**Solution**:
1. Open [fetch-genome.md](fetch-genome.md)
2. Go to **Performance Considerations** → **Parallelization**
3. See recommended concurrency (5-10 parallel downloads)
4. Review **Resource Usage** section

### Use Case 4: Writing Tests

**Scenario**: "How do I test DB_METADATA with test data?"

**Solution**:
1. Open [db-metadata.md](db-metadata.md)
2. Go to **Testing** section
3. Copy unit test code
4. Use provided test data CSV
5. Run and validate results

## 🔍 Search Guide

### Finding Information

#### By Topic

```bash
# Search all module docs for specific topics
cd docs/pipelines/statistics/modules/

# Find BUSCO-related info
grep -r "BUSCO" .

# Find parameter references
grep -r "params\." .

# Find error handling
grep -r "Error" .
```

#### By Module

```bash
# Find all references to a specific module
grep -r "DB_METADATA" .

# Find module outputs
grep -r "output:" .

# Find resource requirements
grep -r "CPUs:" .
```

#### By File Type

```bash
# Find YAML examples
grep -r "```yaml" .

# Find Groovy code examples
grep -r "```groovy" .

# Find bash scripts
grep -r "```bash" .
```

## 💡 Tips & Tricks

### Reading Efficiently

1. **Start with Overview**: Understand what the module does
2. **Check Examples First**: See real usage before diving into details
3. **Use Error Handling**: When troubleshooting, jump straight here
4. **Bookmark Common Modules**: Keep frequently-used docs handy

### Navigation Shortcuts

- **Module Index**: Your navigation hub - use it often
- **Cross-References**: Click "Related Documentation" links
- **Diagrams**: Visual learners - focus on Mermaid diagrams
- **Code Blocks**: Copy-paste examples to learn faster

### Common Patterns

All modules follow similar patterns:

1. **Metadata Propagation**: All modules receive and emit `meta` maps
2. **Version Tracking**: All modules emit `versions.yml`
3. **Error Handling**: Similar error patterns across modules
4. **Resource Labels**: Consistent labeling (python, busco, etc.)

## 📊 Documentation Statistics

### Current Coverage

- **Total Modules**: 14
- **Fully Documented**: 3 (21%)
- **In Progress**: 0
- **Planned**: 11 (79%)

### Documentation Quality

- **Average Length**: 633 lines per module
- **Sections per Module**: ~16
- **Code Examples**: ~10 per module
- **Diagrams**: 1-3 per module

### Content Breakdown

| Section Type | Percentage |
|--------------|------------|
| Examples & Code | 30% |
| Error Handling | 20% |
| Implementation Details | 20% |
| Configuration & I/O | 15% |
| Testing & Troubleshooting | 10% |
| Integration & References | 5% |

## 🛠️ For Contributors

### Documentation Standards

When adding new module documentation:

1. **Follow the Template**: Use existing modules as templates
2. **Include All Sections**: Cover all 16 standard sections
3. **Add Real Examples**: Use actual data, not placeholders
4. **Test Examples**: Ensure all code examples work
5. **Link Related Docs**: Add cross-references
6. **Update Index**: Add module to index.md

### Quality Checklist

- [ ] Overview is clear and concise
- [ ] All inputs/outputs documented
- [ ] Complete metadata map provided
- [ ] Usage example is runnable
- [ ] Error handling covers common cases
- [ ] Performance metrics included
- [ ] Integration diagram added
- [ ] Testing section complete
- [ ] All links work
- [ ] Code blocks have language tags

## 📞 Getting Help

### Where to Look

1. **Module Documentation**: Start with specific module docs
2. **Module Index**: Overview and navigation
3. **Error Handling Sections**: Specific error solutions
4. **Troubleshooting Guides**: Step-by-step debugging

### What's Not Here

This documentation focuses on **individual modules**. For other topics, see:

- **Workflows**: How modules work together (coming soon)
- **Parameters**: Global parameter reference (coming soon)
- **Configuration**: Pipeline configuration (coming soon)
- **Installation**: Setup and requirements (see main README)

## 🎓 Learning Resources

### For Beginners

1. Read [Module Index](index.md) for overview
2. Study [DB_METADATA](db-metadata.md) - simplest module
3. Try [FETCH_GENOME](fetch-genome.md) - file handling
4. Practice with [CLEANING](cleaning.md) - utility patterns

### For Intermediate Users

1. Understand metadata propagation between modules
2. Study channel input/output patterns
3. Learn resource allocation strategies
4. Practice error handling

### For Advanced Users

1. Extend modules with custom features
2. Optimize resource allocation
3. Implement caching strategies
4. Create new modules following patterns

## 🔗 External Resources

### Nextflow

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)
- [nf-core Best Practices](https://nf-co.re/docs/contributing/guidelines)

### Tools

- [BUSCO User Guide](https://busco.ezlab.org/)
- [OMArk Documentation](https://github.com/DessimozLab/OMArk)
- [Ensembl API](https://www.ensembl.org/info/docs/api/)

### Databases

- [Ensembl Database Schema](https://www.ensembl.org/info/docs/api/core/core_schema.html)
- [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly)
- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

## 📝 Quick Reference

### Module Cheat Sheet

| Module | Primary Input | Primary Output | Typical Runtime |
|--------|---------------|----------------|-----------------|
| DB_METADATA | dbname | metadata.txt | 5-15 sec |
| FETCH_GENOME | gca | genome.fna | 5-30 min |
| CLEANING | taxon_id, run_id | (none) | 1-60 sec |

### Common Parameters

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `params.outdir` | Output directory | `./results` |
| `params.cacheDir` | Cache directory | `./cache` |
| `params.mysqlUrl` | Database URL | Required |
| `params.mysqluser` | Database user | Required |

### Resource Labels

| Label | CPUs | Memory | Typical Use |
|-------|------|--------|-------------|
| default | 1 | 2 GB | Light processes |
| python | 2 | 4 GB | Python scripts |
| fetch_file | 2 | 2 GB | File downloads |
| busco | 8 | 16 GB | BUSCO analysis |
| omamer | 4 | 8 GB | OMArk processes |

## 🗺️ Documentation Roadmap

### Next Priorities

1. **BUSCO_DATASET** - Critical for BUSCO workflow
2. **BUSCO_GENOME_LINEAGE** - Core BUSCO functionality
3. **BUSCO_PROTEIN_LINEAGE** - Parallel to genome mode
4. **FETCH_PROTEINS** - Completes data retrieval set

### Future Plans

- Workflow documentation
- Parameter reference guide
- Configuration guide
- FAQ section
- Video tutorials
- Interactive examples

---

## 📍 You Are Here

```
ensembl-genes-nf/
└── docs/
    └── pipelines/
        └── statistics/
            └── modules/
                ├── README.md ← YOU ARE HERE
                ├── index.md
                ├── db-metadata.md
                ├── fetch-genome.md
                └── cleaning.md
```

**Next Steps**: 
1. Read [index.md](index.md) for module overview
2. Explore individual module documentation
3. Try examples with your data
4. Refer back as needed

---

**Last Updated**: 2026-02-06  
**Documentation Version**: 1.0.0  
**Maintained By**: Ensembl Genes Team
