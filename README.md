# crt-tpc-match
Matches SBND CRT tracks to wire hits in order to determine the ratio of real CRT hits to false CRT hits

The main ROOT macro is `wires.cpp`. In ROOT, run `wires(n)` where n is the event number to look at, or -1 to look at all
 events. This macro assumes that there is a file named `hitdumper_tree.root` in the directory it is run from, and that 
 the directory above that contains `WireDumpSBND.txt` and `StripDumpSBND.txt`. Config variables can be found near the 
 top of the file. 
 
 `DRAW`: controls whether the hits, track, etc are drawn.

 `FLIP`: ROOT draws the z-axis as the vertical, setting this to true makes the y-axis the vertical.
 
 `CHITS`, `WHITS`: The max number of crt and wire hits that can be read.
 
### Algorithm
Reconstructs all possible CRT tracks across the top CRT plane(s) to the bottom plane. Determines location of wire hits 
 from intersections of wires with hits that line up in time. Then, for each CRTTrack, "walks" down the track one cm at a
 time and counts the number of nearby wire hits. This score is used to decide which tracks are real and which are noise:
 a track is only marked as real if its score is above a certain threshold (currently 0) and it doesn't share any CRT 
 hits with a higher scored track.