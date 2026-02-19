# RELEASING.md — JaFaR-Suite (GitHub → Zenodo DOI)

This document describes how to publish a new **GitHub Release** of **JaFaR-Suite** and obtain a **citable Zenodo DOI** for that exact version.

Repository: https://github.com/CTP-Lab/JaFaR-Suite  
Zenodo DOI (v1.0.0): https://doi.org/10.5281/zenodo.18695152

---

## 1) One-time setup (only once)

### 1.1 Requirements
- The GitHub repository must be **public** to use the Zenodo GitHub integration.
- Make sure the repository contains (recommended):
  - `README.md`
  - `LICENSE`
  - `CITATION.cff` (and optionally `CITATION.bib`)
  - `RELEASING.md` (this file)

### 1.2 Enable Zenodo GitHub integration
1. Log into Zenodo.
2. Go to **Profile → GitHub**.
3. Click **Sync** (or **Sync now**).
4. Find the repository **CTP-Lab/JaFaR-Suite**.
5. Toggle the switch **ON** to enable archiving.

✅ After this, every new GitHub Release will be automatically archived by Zenodo.

---

## 2) Versioning policy (recommended)

Use Semantic Versioning:

- **MAJOR** `X.0.0` — breaking changes
- **MINOR** `0.X.0` — new features, backward compatible
- **PATCH** `0.0.X` — bug fixes, backward compatible

Examples:
- `v1.0.1` (patch)
- `v1.1.0` (minor)
- `v2.0.0` (major)

---

## 3) Before releasing (pre-release checklist)

### 3.1 Run a quick smoke test
- Launch each app locally to confirm it starts without error:
  - `apps/JaFaR_MZmine_to_W4M_bidirectional`
  - `apps/JaFaR_Untargeted_Metabolomics_Univariate_Pipeline`
  - `apps/JaFaR_Untargeted_Metabolomics_Multivariate_Pipeline`
  - `apps/JaFaR_CSV_and_MGF_filter_just_before_annotation`

### 3.2 Update version metadata (repo)
Before creating the GitHub release, update:
- `CITATION.cff`
  - `version: "X.Y.Z"`
  - `date-released: "YYYY-MM-DD"` (today’s date)
  - (optional) keep `doi:` pointing to the latest published DOI or update it after Zenodo mints the new DOI
- (optional) `CITATION.bib`
  - update `version = {vX.Y.Z}`

> Note: Zenodo will generate a **new DOI per release**. You will update DOI values after Zenodo mints them (Section 5).

### 3.3 Commit everything
Make sure all changes are pushed to GitHub before creating the release.

---

## 4) Create a GitHub Release (this triggers Zenodo)

1. Open: https://github.com/CTP-Lab/JaFaR-Suite
2. Go to **Releases**.
3. Click **Draft a new release**.
4. **Choose a tag**:
   - Enter `vX.Y.Z` (example: `v1.0.1`)
   - Click **Create new tag: vX.Y.Z on publish** (if it doesn’t exist)
5. **Target**: `main`
6. **Release title**: `vX.Y.Z`
7. **Describe this release** (release notes). Recommended template:

   ```text
   ## Changes
   - Added:
   - Fixed:
   - Changed:

   ## Notes
   - (optional) Any migration steps / known issues
   ```

8. Ensure **Set as a pre-release** is **unchecked** (unless it’s explicitly a beta).
9. Click **Publish release**.

✅ Zenodo will now archive this release and mint a DOI.

---

## 5) After releasing (Zenodo DOI + metadata update)

### 5.1 Find the Zenodo record and DOI
1. Go to Zenodo → **Profile → GitHub**
2. Find the new release record for **JaFaR Suite**
3. Open the Zenodo record page
4. Copy the **version DOI** (example format: `10.5281/zenodo.XXXXXXX`)

Zenodo also provides a **Concept DOI** (points to the latest version).  
**Best practice for papers:** cite the **version DOI** used for the analysis.

### 5.2 Update repository citation files with the new DOI
After Zenodo mints the DOI, update:

- `README.md`
  - DOI badge (optional but recommended)
  - “How to cite” section with the new DOI

- `CITATION.cff`
  - update `doi: "10.5281/zenodo.XXXXXXX"`
  - ensure `version: "X.Y.Z"` matches the release

- (optional) `CITATION.bib`
  - update `doi` and `url`

Then commit the changes (this commit is **after** the release; it’s normal).

---

## 6) Recommended citation format

Use the Zenodo **version DOI**:

Fauquet, J. (2026). *JaFaR-Suite* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.18695152

For future releases, replace `v1.0.0` and the DOI accordingly.

---

## 7) Troubleshooting

### 7.1 Zenodo did not archive the release
- Check Zenodo → **Profile → GitHub**:
  - repo is toggled **ON**
  - click **Sync now**
- Confirm the GitHub release is **published** (not draft)
- Wait a few minutes and refresh

### 7.2 Repo does not appear in Zenodo
- Ensure the GitHub repo is **public**
- Click **Sync now** on Zenodo GitHub page
- If the repo is under an organization, ensure permissions allow Zenodo integration

### 7.3 Wrong version/tag format
- Recommended tag format: `vX.Y.Z` (e.g., `v1.2.0`)
- If you created a wrong tag, create a new correct release tag for the next version.

---

## 8) Quick checklist (copy/paste)

Before release:
- [ ] Run smoke tests on the 4 apps
- [ ] Update `CITATION.cff` version + date
- [ ] Commit + push changes

Release:
- [ ] Create GitHub Release with tag `vX.Y.Z`
- [ ] Add release notes
- [ ] Publish release

After release:
- [ ] Get Zenodo version DOI
- [ ] Update README (badge + “How to cite”)
- [ ] Update `CITATION.cff` (`doi:`)
- [ ] (Optional) Update `CITATION.bib`
- [ ] Commit + push DOI updates
