# Nextflow Linting Guide

This repository uses automated Nextflow linting to maintain code quality and consistency.

## 🔍 Automatic Linting (CI)

### What Gets Checked

The `nextflow-lint.yml` workflow runs automatically on:
- Every push to any branch
- Every pull request

It checks all `.nf` files for:
- Syntax errors
- Code style issues
- DSL2 compliance
- Strict syntax compatibility

## 🛠️ Available Workflows

### 1. `nextflow-lint.yml` - Continuous Linting
- **Trigger**: Automatic on push/PR
- **Purpose**: Check code quality
- **Action**: Reports errors but doesn't modify files
- **Exit**: Fails if linting errors are found

### How to Fix Lint Errors

If the CI check fails, you can auto-format locally:

```bash
# Install or update to Nextflow 25.10.0
export NXF_VER=25.10.0
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Format all files
nextflow lint -format .

# Or format with custom spacing
nextflow lint -format -spaces 4 .

# Commit the changes
git add -A
git commit -m "style: format Nextflow files"
git push
```

## 📋 Lint Command Reference

### Check Only (No Changes)
```bash
nextflow lint .                    # Full output
nextflow lint -o concise .         # Concise output
nextflow lint -o json .            # JSON output (for parsing)
nextflow lint main.nf              # Lint specific file
```

### Format Files
```bash
nextflow lint -format .            # Format with 4 spaces (default)
nextflow lint -format -spaces 2 .  # Format with 2 spaces
nextflow lint -format -tabs .      # Format with tabs
```

### Advanced Options
```bash
# Exclude specific patterns
nextflow lint -exclude ".git,.nf-test,work" .

# Sort declarations
nextflow lint -format -sort-declarations .

# Lint specific directory
nextflow lint pipelines/statistics/
```

## 🔧 Version Compatibility

### Nextflow 25.10.0+ (Current)
- Supports section label syntax for workflow handlers:
  ```groovy
  workflow {
      main:
      // code
      
      onComplete:
      log.info("Done!")
  }
  ```

### Nextflow 25.04.x (Legacy)
- Requires assignment syntax:
  ```groovy
  workflow {
      main:
      // code
  }
  
  workflow.onComplete = {
      log.info("Done!")
  }
  ```

### Downgrade Formatting (if needed)
If you need backward compatibility with Nextflow < 25.10:

```bash
# Use Nextflow 25.04 for formatting
NXF_VER=25.04.8 nextflow lint -format .
```

## 🚫 Disabling Strict Syntax (Not Recommended)

Strict syntax is **enabled by default** in the workflows. If you absolutely need to disable it:

```bash
# Temporarily disable strict parsing
unset NXF_SYNTAX_PARSER

# Or explicitly use old parser
export NXF_SYNTAX_PARSER=v1
```

**Note**: This is not recommended as strict syntax will become mandatory in future Nextflow versions.

## 📖 Common Issues and Fixes

### Issue: "Invalid workflow definition"
**Cause**: Using section labels (`onComplete:`) with Nextflow < 25.10

**Fix**: Either upgrade to 25.10+ or use assignment syntax:
```groovy
workflow.onComplete = { /* code */ }
```

### Issue: "Implicit closure parameter"
**Cause**: Using implicit `it` parameter

**Fix**: Declare parameters explicitly:
```groovy
// ❌ Bad
ch.map { it.toUpperCase() }

// ✅ Good
ch.map { v -> v.toUpperCase() }
```

### Issue: "channel vs Channel"
**Cause**: Using old `Channel` type instead of `channel` namespace

**Fix**: Use lowercase namespace:
```groovy
// ❌ Bad
Channel.of(1, 2, 3)

// ✅ Good
channel.of(1, 2, 3)
```

### Issue: "Process env must be quoted"
**Cause**: Unquoted environment variable names

**Fix**: Quote the variable name:
```groovy
// ❌ Bad
env FOO

// ✅ Good
env 'FOO'
```

## 🎯 Best Practices

1. **Run lint before committing**:
   ```bash
   nextflow lint .
   ```

2. **Use auto-formatting**:
   ```bash
   nextflow lint -format .
   ```

3. **Check specific files during development**:
   ```bash
   nextflow lint pipelines/statistics/main.nf
   ```

4. **Use concise output for quick checks**:
   ```bash
   nextflow lint -o concise .
   ```

## 🔗 Resources

- [Nextflow Lint Documentation](https://www.nextflow.io/docs/latest/reference/cli.html#lint)
- [Strict Syntax Guide](https://www.nextflow.io/docs/latest/strict-syntax.html)
- [DSL2 Migration Guide](https://www.nextflow.io/docs/latest/dsl2.html)
- [Nextflow Style Guide](https://www.nextflow.io/docs/latest/developer/style.html)