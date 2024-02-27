#!/usr/bin/awk -f

{
    if ($0 ~ /^S/) {
        print ">"$2 "\n" $3
    }
}