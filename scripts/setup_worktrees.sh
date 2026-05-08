#!/usr/bin/env bash
# setup_worktrees.sh — create 7 git worktrees for parallel agent work
#
# Run from repo root:
#   bash scripts/setup_worktrees.sh
#
# Result: 7 sibling directories beside the main repo, each on its own
# feature branch, each containing a full checkout. Open one IDE window
# per worktree, dispatch a Claude agent in each.

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
PARENT_DIR="$(dirname "$REPO_ROOT")"
BASE_BRANCH="active-learning-pipeline"

# Sanity: must be on the base branch with clean tree
current_branch="$(git -C "$REPO_ROOT" branch --show-current)"
if [[ "$current_branch" != "$BASE_BRANCH" ]]; then
  echo "ERROR: must run from $BASE_BRANCH (currently on $current_branch)" >&2
  exit 1
fi
if [[ -n "$(git -C "$REPO_ROOT" status --porcelain)" ]]; then
  echo "ERROR: working tree is dirty. Commit or stash first." >&2
  exit 1
fi

# Module → branch name mapping
declare -a MODULES=(
  "00-raw-data"
  "01-data-ground-truth"
  "02-annotation"
  "03-rbp-identification"
  "04-protein-embedding"
  "05-structure-prediction"
  "06-uncertainty-model"
)

echo "Repo root:     $REPO_ROOT"
echo "Parent dir:    $PARENT_DIR"
echo "Base branch:   $BASE_BRANCH"
echo ""

for mod in "${MODULES[@]}"; do
  branch="agent-${mod}"
  worktree_path="${PARENT_DIR}/agent-${mod}"

  if git -C "$REPO_ROOT" worktree list --porcelain | grep -q "^worktree $worktree_path\$"; then
    echo "[SKIP] $worktree_path already exists"
    continue
  fi

  if git -C "$REPO_ROOT" show-ref --quiet "refs/heads/$branch"; then
    echo "[INFO] Branch $branch exists; checking it out"
    git -C "$REPO_ROOT" worktree add "$worktree_path" "$branch"
  else
    git -C "$REPO_ROOT" worktree add "$worktree_path" -b "$branch"
  fi
  echo "[OK]   $worktree_path  →  $branch"
done

echo ""
echo "Done. Worktrees:"
git -C "$REPO_ROOT" worktree list

cat <<'EOF'

==============================================================
Next steps for Alex:
==============================================================
1. Open 7 separate IDE windows, one per worktree directory.
2. In each window, launch Claude with --dangerously-skip-permissions.
3. The agent's first instruction should be:
     "Read INTERFACE.md, LAGUNA.md, and your module's
      <module>/AGENT_TODO.md, then execute the plan."
4. Agents commit locally to their feature branch. Tomorrow morning,
   review the 7 branches and merge the good ones into
   active-learning-pipeline.

To remove a worktree later (cleanup):
   git worktree remove ../agent-<module>
   git branch -d agent-<module>
EOF
