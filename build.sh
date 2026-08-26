#!/bin/sh

for DOCUMENT in `ls src/*.adoc`; do
	OUTPUT="$(basename $DOCUMENT .adoc).html"
	[ $(stat $DOCUMENT --printf="%X") -lt $(stat $OUTPUT --printf="%X") ] || asciidoctor $DOCUMENT -r asciidoctor-diagram -o $OUTPUT
done
