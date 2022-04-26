#!/bin/bash
cat << EOF
inputRedirect = {};
inputRedirect["01"]="Universal"; // genetic code
inputRedirect["02"]="$(readlink -f $1)"; // codon data
inputRedirect["03"]="$(readlink -f $2)"; // tree
inputRedirect["04"]="${3:-All}"; // Test for selection on a branch
inputRedirect["05"]=""; // complete selection


ExecuteAFile ("/usr/local/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf", inputRedirect);
EOF
