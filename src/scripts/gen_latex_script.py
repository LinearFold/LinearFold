#!/usr/bin/env python2.7

import sys
import math
from collections import defaultdict
logs = sys.stderr

preamble = r'''
\documentclass{standalone}
%\usepackage{fullpage}
\pagestyle{empty}
\usepackage{tikz}
\usepackage{tkz-euclide}
\usepackage{siunitx}
\usetikzlibrary{shapes, shapes.multipart}

\usepackage{verbatim}
\usepackage{lipsum}

\begin{document}
'''

picturepre = r'''
%\hspace{-3cm}
%\resizebox{1.2\textwidth}{!}{
\begin{tikzpicture}[darkstyle/.style={}, scale=2] %circle,draw,fill=gray!10}]
'''

print preamble

dataset = "."
MAXLEN = 5650
MINLEN = 50
circular = True

lbs = [
    '(',
    '[',
    '{',
    '<'
]

rbs = [
    ')',
    ']',
    '}',
    '>'
]

counter_clockwise = False # hzhang
rotate = 180 - 4.5


def drawarc_clockwise(a, b, deg, style, length, lengthfix):
    # a,b = b,a
    angle_a = 360./(length+lengthfix)*a
    angle_b = 360./(length+lengthfix)*b
    angle_a = angle_a-70 + rotate
    angle_b = angle_b-70 + rotate

    alpha = (angle_b-angle_a)
    pa = angle_a * math.pi / 180.
    pb = angle_b * math.pi / 180.
    palpha = alpha * math.pi / 180.
    if alpha < 170:
        print "\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_b,
                                                                                                    angle_b + 90,
                                                                                                    angle_a + 270,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(palpha/2.),
                                                                                                    a, b, length
        )
    elif alpha > 190:
        beta = 360. - alpha
        pbeta = beta * math.pi / 180.
        print "\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_a,
                                                                                                    angle_a + 90,
                                                                                                    angle_b - 90,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(pbeta/2.),
                                                                                                    a, b, length
        )
    else:
        print "\\draw[line width = 0.0mm] %s (%d) to [bend left=%.1f] (%d);" % (style, 
                                                                              a, 
                                                                              2*deg,
                                                                              b)


def drawarc_counterclockwise(a, b, deg, style, length, lengthfix): ## counterclock wise
    a,b = b,a # hzhang
    angle_a = 360./(length+lengthfix)*a
    angle_b = 360./(length+lengthfix)*b
    alpha = (angle_b-angle_a)
    pa = angle_a * math.pi / 180.
    pb = angle_b * math.pi / 180.
    palpha = alpha * math.pi / 180.
    if alpha < 170:
        print "\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_b,
                                                                                                    angle_b + 90,
                                                                                                    angle_a + 270,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(palpha/2.),
                                                                                                    a, b, length
        )
    elif alpha > 190:
        beta = 360. - alpha
        pbeta = beta * math.pi / 180.
        print "\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_a,
                                                                                                    angle_a + 90,
                                                                                                    angle_b - 90,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(pbeta/2.),
                                                                                                    a, b, length
        )
    else:
        print "\\draw[line width = 0.0mm] %s (%d) to [bend left=%.1f] (%d);" % (style, 
                                                                              a, 
                                                                              2*deg,
                                                                              b)


def agree(pres, pref, a, b): ## pres[a] = b
    if pref[a] == b:
        return True
    elif pref.get(a-1,-1) == b or pref.get(a+1,-1) == b:
        return True
    elif pref.get(b-1,-1) == a or pref.get(b+1,-1) == a:
        return True
    else:
        return False

note = False
notes = ""
num_hasknot = 0
for index, line in enumerate(sys.stdin):
    pairs = []
    goldpairs = []
    pairset = set()
    respair = defaultdict(lambda: -1)
    goldpair = defaultdict(lambda: -1)
    bases = ['']

    pseudoknot = 0
    tmp = line.strip().split()
    if len(tmp) == 3:
        seq, res, ref = tmp
        filename = "seq %d" % index
    elif len(tmp) == 4:
        note = True
        seq, res, ref, notes = tmp
    elif len(tmp) == 5:
        note = True
        seq, res, ref, filename, notes = tmp
    else:
        print "input format error!"
        sys.exit(1)

    stacks = []
    for _ in xrange(len(lbs)):
        stacks.append([])
    for i, item in enumerate(res):
        if item in lbs:
            stackindex = lbs.index(item)
            stacks[stackindex].append(i)
        elif item in rbs:
            stackindex = rbs.index(item)
            left = stacks[stackindex][-1]
            stacks[stackindex] = stacks[stackindex][:-1]
            pairs.append((left+1,i+1, stackindex))
            respair[left+1] = i+1
            respair[i+1] = left+1
            pairset.add((left+1,i+1))
    notes += ";pair=%d" % (len(respair)//2)

    stacks = []
    for _ in xrange(len(lbs)):
        stacks.append([])
    for i, item in enumerate(ref):
        if item in lbs:
            stackindex = lbs.index(item)
            stacks[stackindex].append(i)
        elif item in rbs:
            stackindex = rbs.index(item)
            left = stacks[stackindex][-1]
            stacks[stackindex] = stacks[stackindex][:-1]
            goldpairs.append((left+1,i+1, stackindex))
            goldpair[left+1] = i+1
            goldpair[i+1] = left+1

    length = len(seq)
    bases = [''] + list(seq)


    if length > MAXLEN:
        print >> logs, "too long (%d %s). Sequence length should in range [50, 5650]" % (length, "nt")
        continue # too long    

    if length < MINLEN:
        print >> logs, "too short (%d %s). Sequence length should in range [50, 5650]" % (length, "nt")
        continue # too short    

    print picturepre


    lengthfix = int(length/9.0)
    for i, base in enumerate(bases[1:], 1):
        if circular:
            angle = 360./(length+lengthfix)*(i)
            
            if counter_clockwise: # hzhang
                angle = 450-angle -20
            else:
                angle = angle-70 + rotate 
            
            print "\\node [darkstyle]  (%d) at (%f:10cm) {};" % (i, angle)
            # print len(seq)
            if length <= 100:
                gap = 5
            elif length <= 200:
                gap = 10
            elif length <= 300:
                gap = 20
            elif length <= 700:
                gap = 50
            elif length <= 2000:
                gap = 100
            elif length <= 3000:
                gap = 200
            elif length <= 5000:
                gap = 300
            else:
                gap = 400

            if i % gap == 0 and i < len(bases)-10:
                print "\\node [scale=2]           (%d,1) at (%f:10.8cm) {\Huge %d};" % (i, angle, i)
            if i > 1:
                print "\\draw (%d.center) -- (%d.center);" % (i, i-1)
      
        else:
            print "\\node [darkstyle]  (%d) at (%d,0) {};" % (i, i)
            if i % 5 == 0:
                print "\\node []           (%d,1) at (%d,-1) {%d};" % (i, i, i)

    if circular:

        for j in range(4):
            i += 1
            angle = 360./(length+lengthfix)*(i)
            angle = 450-angle-20
            print "\\node [darkstyle]  (%d) at (%f:10cm) {};" % (i, angle)
            # print "\\draw (%d.center) -- (%d.center);" % (i, i-1)

        angle = 360./(length+lengthfix)
        angle = 450-angle
        if length <= 100:
            angle5 = angle - 20
            angle3 = angle + 20
        if length <= 200:
            angle5 = angle - 17
            angle3 = angle + 17
        else:
            angle5 = angle - 13
            angle3 = angle + 13
        # print "\\node [scale=2](%d,1) at (%f:10cm) {\LARGE \\textbf{%s}};" % (i, angle5, "5'")
        if counter_clockwise:
            print "\\node [scale=2](3prime) at (%f:10cm) {\LARGE \\textbf{%s}};" % (angle3, "3'")
            print "\\node [right=9.5cm of 3prime, scale=2] {\LARGE \\textbf{%s}};" % "5'"
        else:
            print "\\node [scale=2](3prime) at (%f:10cm) {\LARGE \\textbf{%s}};" % (angle3, "5'")
            print "\\node [right=9.5cm of 3prime, scale=2] {\LARGE \\textbf{%s}};" % "3'"

    for a, b, stackindex in goldpairs:
        if agree(goldpair, respair, a, b):
            ifdraw = False
        else:
            ifdraw = True

        if not ifdraw:
            continue
        color = "gray!20"

        pair = (bases[a]+bases[b]).upper()
        style = "[" + color

        if stackindex == 0:
            if b-a < length+lengthfix-b+a:
                style += ",thick"
            else:
                # style += ",dashed"
                style += ",thick"
        else:
            # style += ",dashed"
            style += ",thick"
        style += "]"

        if circular: 
            dist = b - a
            revdist = length+lengthfix - dist

            deg = 90 * (0.5-(dist+.0) / (length+lengthfix+.0))

        else:
            deg = 20 # linear layout

        if counter_clockwise:
            drawarc_counterclockwise(a, b, deg, style, length, lengthfix)
        else:
            drawarc_clockwise(a, b, deg, style, length, lengthfix)

    for a, b, stackindex in pairs:
        if stackindex > 0:
            color = "red"
        else:
            if agree(respair, goldpair, a, b):
                color = "blue"
            else:
                color = "red"

        pair = (bases[a]+bases[b]).upper()
        style = "[" + color

        if stackindex == 0:
            if b-a < length+lengthfix-b+a:
                style += ",thick"
            else:
                # style += ",dashed"
                style += ",thick"
        else:
            # style += ",dashed"
            style += ",thick"
        style += "]"

        if circular: 
            dist = b - a
            revdist = length+lengthfix - dist

            deg = 90 * (0.5-(dist+.0) / (length+lengthfix+.0))

        else:
            deg = 20 # linear layout

        if counter_clockwise:
            drawarc_counterclockwise(a, b, deg, style, length, lengthfix)
        else:
            drawarc_clockwise(a, b, deg, style, length, lengthfix)

    # if note:
    #     print "\\node[align=center,font=\\bfseries, yshift=2em] (title) at (current bounding box.north) {\Huge %s};" % (notes)
    print "\\end{tikzpicture}"
    # print r"} \newpage " # resizebox

print "\\end{document}"

print >> logs, "%d out of %d sequences have pseudoknots" % (num_hasknot, index) #TODO
