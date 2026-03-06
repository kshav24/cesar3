# CESAR3 Agent Rules

## CRITICAL INSTRUCTION 1
Never run cat inside a bash command to create or append to files.
Never use ls for listing, grep inside bash, or sed for replacing.
Use the provided file tools instead.

## CRITICAL INSTRUCTION 2
One step. One commit. One diff. Then STOP and wait.
After every commit run: git diff HEAD~1
Output the full diff and wait for explicit "approved" or "next step"
before proceeding. Never self-approve.

## CRITICAL INSTRUCTION 3
After every commit, update docs/changes.tex in the SAME commit.
Add a new \section{} with:
- \subsection{What Changed}: file name + exact diff hunks
- \subsection{Why}: the architectural reason, not just what
- \subsection{Code}: \begin{lstlisting}[language=C] block
- \subsection{Status}: PLANNED / IN PROGRESS / COMMITTED / TESTED
Never skip the LaTeX update. Always same commit as the code change.
