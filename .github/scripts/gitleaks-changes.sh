#!/usr/bin/env bash
# Scan only the changes being introduced for secrets, not content that was already
# scanned elsewhere. Run by the gitleaks hook in .pre-commit-config.yaml.
#
# - GITLEAKS_BASE set (CI): the commits in GITLEAKS_BASE..HEAD, plus the lines that
#   merge commits in that range introduced themselves (e.g. conflict resolutions),
#   not what they brought in from the merged branch. Set but empty (a manual
#   workflow_dispatch run has no base): the whole history.
# - Merge commit in progress: what differs from the branch being merged in (this
#   branch's changes and the conflict resolutions), not everything merged in.
# - Otherwise: the staged changes (the upstream hook's behavior).
set -euo pipefail

# Lines of a combined diff (git show --cc) that are in none of the merge's parents,
# with the file headers. With N parents those lines start with N "+".
merge_introduced_lines() {
    local merge=$1 parents
    parents=$(($(git rev-list --parents -n 1 "$merge" | wc -w) - 1))
    git show --cc --format= "$merge" | awk -v n="$parents" '
        BEGIN { plus = sprintf("%" n "s", ""); gsub(/ /, "+", plus) }
        /^diff --cc / { print; next }
        /^\+\+\+ (b\/|\/dev\/null)/ { next }
        substr($0, 1, n) == plus { print }
    '
}

if [[ -n "${GITLEAKS_BASE+set}" ]]; then
    if [[ -z "${GITLEAKS_BASE}" ]]; then
        echo "gitleaks: no base given, scanning the whole history"
        exec gitleaks git --redact --verbose
    elif git rev-parse --quiet --verify "${GITLEAKS_BASE}^{commit}" >/dev/null; then
        range=("${GITLEAKS_BASE}..HEAD")
    else
        # e.g. first push of a branch (all zeros) or a base that was force-pushed away
        echo "gitleaks: base ${GITLEAKS_BASE} not found, scanning the last commit only" >&2
        range=(-1 HEAD)
    fi

    status=0
    echo "gitleaks: scanning commits ${range[*]}"
    gitleaks git --redact --verbose --log-opts="${range[*]}" || status=1

    # git log -p (used by gitleaks git) shows no diff for merge commits
    for merge in $(git rev-list --merges "${range[@]}"); do
        echo "gitleaks: scanning the lines merge commit ${merge} introduced"
        merge_introduced_lines "$merge" | gitleaks stdin --redact --verbose || status=1
    done
    exit "$status"
fi

if git rev-parse --quiet --verify MERGE_HEAD >/dev/null; then
    # Added lines only (with the file headers), removed lines are not introduced
    git diff --cached --unified=0 MERGE_HEAD | { grep -E '^(diff --git |\+)' || true; } | gitleaks stdin --redact --verbose
    exit
fi

exec gitleaks git --pre-commit --redact --staged --verbose
