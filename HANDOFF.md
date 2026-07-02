# InstrumentQC — session handoff (2026-07-02)

Context for resuming work in a new Claude Code session. (This file is a scratch
note — delete it whenever. It is not in the Quarto render whitelist.)

Durable project facts (layout, render commands, gotchas) now live in
[CLAUDE.md](CLAUDE.md) — read that first, it auto-loads for any Claude Code
session in this repo. Cross-session diagnostic history (known-bad runs, fit
caveats, methodology provenance) lives in this account's memory files — see
`instrumentqc-project.md`, `instrumentqc-render-gotchas.md`,
`evo-8peak-known-runs-and-qb-caveats.md`, `evo-6peak-urcp-pipeline.md`,
`qb-sensitivity-reference.md`.

---

## Current state (as of this session)
Branch `main`, otherwise clean. This session:
- Added `CLAUDE.md` (consolidated durable project facts out of memory).
- Added `.claude/settings.json` (read-only permission allowlist for
  `preview_screenshot`/`preview_list` — see `fewer-permission-prompts` skill
  output; most common commands were already covered by Claude Code's built-in
  auto-allow list, so this only added 2 entries).
- Staged deletion of `pages/AuroraEVO_8Peak_previous.qmd` (comparison artifact,
  no longer needed).

Nothing else in flight from a prior session.

## Next steps / open items
1. Pre-existing loose ends (untouched, still pending the user's call):
   `README.md` modified, `My Working Index.qmd` deleted. Help-page screenshots
   are still the original Maryland ones (could be replaced with screenshots of
   this dashboard).
2. Unanswered question from the user: whether their claude.ai-account
   "Instructions for Claude" apply in Claude Code (they generally do NOT — Code
   uses CLAUDE.md + settings, not the web-app personal instructions).
