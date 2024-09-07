#!/bin/sh

for DOCUMENT in `ls src/*.adoc`; do
	OUTPUT="$(basename $DOCUMENT .adoc).html"
	if [ $(stat $DOCUMENT --printf="%X") -gt $(stat $OUTPUT --printf="%X") ]; then
		asciidoctor $DOCUMENT -o $OUTPUT
	fi
done
