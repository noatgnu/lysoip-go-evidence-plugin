# Examples

Sample input files for manual testing and for CI.

`params.json` holds the `--param`/`--params-json` values `cauldron job run` uses against this
plugin's inputs (see `plugin.yaml`'s `inputs:` section, or run
`cauldron plugin inputs lysoip-go-evidence` once installed). CI runs this automatically via
`.github/workflows/test-plugin.yml`.

```json
{
  "gene_lookup_file": "examples/gene_lookup.tsv"
}
```
