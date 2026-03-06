---
description: Updates the CESAR3 change log LaTeX file after every code commit. Use this skill whenever a C file is modified and committed in the cesar3 repo.
---

# CESAR LaTeX Change Log Updater

After every code commit to the cesar3 repository, append a new entry to docs/changes.tex following this exact template:

\section{Step N: <short title>}
\subsection{What Changed}
File \texttt{<filename>} was modified. The diff is:
\begin{lstlisting}[language=C]
<paste exact diff hunk here>
\end{lstlisting}

\subsection{Why}
<architectural reason — must reference the parallel rewrite goals>

\subsection{Status}
COMMITTED

Then: git add docs/changes.tex && git add <code file> && git commit --amend --no-edit
This amends the code commit so LaTeX and code are always in the same commit.
