import json
import sys
import math
import subprocess
import os

from drawtetrado.svg_painter import Point, ConnType, ConnFlow

class Nucleotide:
    def FindConnections(self, used):
        # Find to what it is connected.
        connected_to = ""
        curr_min = sys.maxsize
        for name, nucl in used.items():
            if nucl["index"] > self.index and \
               curr_min > nucl["index"] and \
               self.chain == nucl["chain"]:
                curr_min = nucl["index"]
                connected_to = name
        return connected_to

    def Block(self, width, height, angle):
        self.coords = []
        self.coords.append(Point(0.0, 0.0))
        if math.tan(math.radians(angle)) != 0:
            shift = height / math.tan(math.radians(angle))
        else:
            shift = height / math.tan(math.radians(0.01))
        self.coords.append(Point(shift, -height))
        self.coords.append(Point(shift + width, -height))
        self.coords.append(Point(width, 0.0))
        self.center = Point((shift + width) / 2.0, -height / 2.0)

    def ShiftBlock(self, shift):
        self.coords[0] += shift
        self.coords[1] += shift
        self.coords[2] += shift
        self.coords[3] += shift
        self.center += shift

    def CalculateCoordinates(self, conf):
        shift = Point(0, 0)
        tan_val = math.tan(math.radians(conf.angle))
        sin_val = math.sin(math.radians(conf.angle))
        cos_val = math.cos(math.radians(conf.angle))
        tan_val_rev = math.tan(math.radians(90.0 - conf.angle))

        shift.y -= self.tetrade_no * (sin_val * (conf.longer + conf.shorter + conf.spacing) + \
                                      conf.tetrade_spacing)
        if self.GetOnzPlusMinus() == "-":
            width_0 = conf.longer
            height_0 = sin_val * conf.shorter
            width_1 = conf.shorter
            height_1 = sin_val * conf.longer
        else: # Should be just this: self.onz == "?+": 
            width_0 = conf.shorter
            height_0 = sin_val * conf.longer
            width_1 = conf.longer
            height_1 = sin_val * conf.shorter

        if self.position == 0:
            self.Block(width_0, height_0, conf.angle)
        elif self.position == 1:
            self.Block(width_1, height_1, conf.angle)
            shift.y -= (width_1 + conf.spacing) * sin_val
            shift.x += (width_1 + conf.spacing) * cos_val
        elif self.position == 2:
            self.Block(width_0, height_0, conf.angle)
            shift.y -= (width_0 + conf.spacing) * sin_val
            shift.x += (width_1 + conf.spacing) + \
                       ((width_0 + conf.spacing) * sin_val) * tan_val_rev
        elif self.position == 3:
            self.Block(width_1, height_1, conf.angle)
            shift.x += width_0 + conf.spacing

        self.ShiftBlock(shift)

    def GetOnzPlusMinus(self):
        if self.onz_full[-1] == "-":
            return "-"
        return "+"

    def GetOnz(self):
        res = self.onz_full[0:1].lower() + "_"
        if self.GetOnzPlusMinus() == "-":
            res += "minus"
        else:
            res += "plus"
        # Only o/n/z_plus/minus are allowed
        if res != "o_plus" and res != "o_minus" and \
           res != "n_plus" and res != "n_minus" and \
           res != "z_plus" and res != "z_minus":
            return "onz_default"
        return res

    def __init__(self, data, used_nucl, tetr_no, tetr_onz, pos):
        self.number = data["number"]
        self.short_name = data["shortName"]
        self.full_name = data["fullName"]
        self.chain = data["chain"]
        self.index = data["index"]
        self.onz_full = tetr_onz
        self.onz = self.GetOnz()
        self.bond = data["glycosidicBond"]
        self.tetrade_no = tetr_no
        self.position = pos
        self.connected_to = self.FindConnections(used_nucl)
        self.connected_from = ""
        self.coords = []
        self.center = Point(0, 0)
        self.connection_type = ConnType.UNKNOWN

        # Loop
        self.connection_label = 0
        self.connection_label_type = ""

        # Parameters regarding connection from this nucleotide to another.
        # out - flow starting the connection
        # in - flow ending the conneciton
        self.flow_out = ConnFlow.UNKNOWN
        self.flow_in = ConnFlow.UNKNOWN


        # Priority of drawing.
        self.priority_conn = -1
        self.priority_edge = -1
        self.priority_nucl = -1

    def UpdatePriorities(self, nucleotides):
        # 1
        if self.connection_type == ConnType.LEFT or \
           self.connection_type == ConnType.LEFT_CROSS or \
           self.connection_type == ConnType.FRONT_TO_BACK:
            self.priority_conn = 1
        elif self.connection_type == ConnType.FRONT_BACK_CROSS and \
             (self.position == 1 or self.position == 2):
            self.priority_conn = 1
        elif self.connection_type == ConnType.SAME_LEVEL and \
             (self.flow_out == ConnFlow.DOWN and self.flow_in == ConnFlow.DOWN):
            self.priority_conn = 1
        elif self.connection_type == ConnType.SIMPLE and self.position == 1:
            self.priority_conn = 1
        # 2 + 3
        if self.connection_type == ConnType.SAME_LEVEL and \
             (self.flow_out == ConnFlow.UP and self.flow_in == ConnFlow.DOWN):
            self.priority_nucl = 2
            self.priority_conn = 3
        elif self.connection_type == ConnType.SAME_LEVEL and \
             (self.flow_out == ConnFlow.DOWN and self.flow_in == ConnFlow.UP):
            if self.connected_to != "":
                nucleotides[self.connected_to].priority_nucl = 2
            self.priority_conn = 3
        # 4
        if self.priority_nucl == -1:
            self.priority_nucl = 4
        # 5 (Have to be later to be on top of the nucleotides)
        if self.priority_edge == -1:
            self.priority_edge = 5
        # 6 (all other connections are processed)
        if  self.priority_conn == -1:
            self.priority_conn = 6


class Quadruplex:

    
    def AddLoopBulge(self, structure):

        nucl_names = set(self.nucl_quad.keys())

        for quadruplex_loops in self.loops:

            for loop in quadruplex_loops:
                
                loop_nucleotides = loop.get("nucleotides", [])
                if len(loop_nucleotides) == 0:
                    continue


                loop_start = structure.nucleotides.get(loop_nucleotides[0])
                loop_finish = structure.nucleotides.get(loop_nucleotides[-1])

                prev_name = ""
                next_name = ""
                prev_index = -1
                next_index = sys.maxsize

                for name, nucl in self.nucl_quad.items():
                    if nucl.chain != loop_start["chain"]:
                        continue

                    if nucl.index < loop_start["index"] and nucl.index > prev_index:
                        prev_index = nucl.index
                        prev_name = name

                    if nucl.index > loop_finish["index"] and nucl.index < next_index:
                        next_index = nucl.index
                        next_name = name


                if prev_name != "" and next_name != "":
                    if self.nucl_quad[prev_name].connected_to == next_name:
                        self.nucl_quad[prev_name].connection_label = len(loop_nucleotides)
                        self.nucl_quad[prev_name].connection_label_type = "loop"


        # BULGE
        for tetrad_tracts in self.tracts:
            for tract in tetrad_tracts:
                for i in range(len(tract) - 1):
                    a = tract[i]
                    b = tract[i + 1]



                    nucl_a = self.nucl_quad[a]
                    nucl_b = self.nucl_quad[b]

                    if nucl_a.chain != nucl_b.chain:
                        continue


                    id_start = min(nucl_a.index, nucl_b.index)
                    id_finish = max(nucl_a.index, nucl_b.index)
                    bulge_count = 0


                    for name, data in structure.nucleotides.items():
                        if data["chain"] == nucl_a.chain and \
                        data["index"] > id_start and \
                        data["index"] < id_finish and \
                        name not in nucl_names:
                            bulge_count += 1

                    if bulge_count == 0:
                        continue

                    if nucl_a.connected_to == b:
                        label_nucl = nucl_a
                    elif nucl_b.connected_to == a:
                        label_nucl = nucl_b
                    else:
                        continue

                    if label_nucl.connection_label == 0:
                        label_nucl.connection_label = bulge_count
                        label_nucl.connection_label_type = "bulge"



    def Rotate(self, shift):
        shift %= 4

        if shift == 0:
            return

        for nucl in self.nucl_quad.values():
            nucl.position = (nucl.position + shift) % 4


    

    def ForceFivePrime(self):
        first_chain = next(iter(self.chains.values()))
        five_name = first_chain["first"]
        five_prime = self.nucl_quad[five_name]

        max_level = len(self.tetrads) - 1

        # 5' na górze
        if five_prime.tetrade_no == max_level:
            for nucl in self.nucl_quad.values():
                nucl.tetrade_no = max_level - nucl.tetrade_no
            self.tetrads.reverse()
            self.tracts.reverse()

            five_prime = self.nucl_quad[five_name]


        # 5' w środku
        if 0 < five_prime.tetrade_no < max_level:
            five_prime = self.nucl_quad[five_name]


        curr_pos = five_prime.position
        shift = (3 - curr_pos) % 4
        self.Rotate(shift)

        new_tetrads = [["", "", "", ""] for _ in range(len(self.tetrads))]
        for name, nucl in self.nucl_quad.items():
            new_tetrads[nucl.tetrade_no][nucl.position] = name
        
        self.tetrads = new_tetrads


        
    
    def UsedNucleotides(self, tetrad, nucl):
        used = {}
        for name, data in tetrad.items():
            used[data["nt1"]] = nucl[data["nt1"]]
            used[data["nt2"]] = nucl[data["nt2"]]
            used[data["nt3"]] = nucl[data["nt3"]]
            used[data["nt4"]] = nucl[data["nt4"]]
        return used

    def PrepareNucleotides(self, structure, quadruplex_id, tetrad_id):
        nucl = structure.nucleotides
    
        if tetrad_id >= 0:
            tetrads = structure.single_tetrads[quadruplex_id][tetrad_id]
        else:
            tetrads = structure.tetrads[quadruplex_id]
    
        tetrads_order = structure.tetrads_order[quadruplex_id]
        used_nucl = self.UsedNucleotides(tetrads, nucl)
        pairs = structure.basePairs

        
        link_map = {}
        relation_map = {}
        
        for p in pairs:
            if not p.get("inTetrad", False):
                continue
                
            lw = p.get("lw", "")

            if lw not in ("cWH", "cHW"):
                continue
                
    
            a = p["nt1"]
            b = p["nt2"]

            link_map.setdefault(a,set()).add(b)
            link_map.setdefault(b,set()).add(a)
            
            relation_map[(a,b)] = lw
            if lw == "cWH":
                relation_map[(b,a)] = "cHW"
            else:
                relation_map[(b,a)] = "cWH"
    
        
        tetr_no = 0
    
        for tetrad_name in tetrads_order:
            if tetrad_name not in tetrads:
                continue
    
            tetrad = tetrads[tetrad_name]
    
            nt1 = tetrad["nt1"]
            nt2 = tetrad["nt2"]
            nt3 = tetrad["nt3"]
            nt4 = tetrad["nt4"]
            onz = tetrad["onz"]
    
            cycle = [nt1, nt2, nt3, nt4]
    

            order = [nt1]
            used = {nt1}
    
            while len(order) < 4:
                last = order[-1]
                found = False
    
                neighbors = link_map.get(last, ())
    
 
                for nt in sorted(neighbors):
                    if nt in cycle and nt not in used:
                        order.append(nt)
                        used.add(nt)
                        found = True
                        break

                if not found:
                    break

            
            if len(order) != 4:
                order = [nt1, nt2, nt3, nt4]

            lw_pair = relation_map.get((order[0], order[1]),"")

            if onz.endswith("+"):
                expected_lw = "cWH"
            else:
                expected_lw = "cHW"

            if lw_pair != "" and lw_pair != expected_lw:
                order = [order[0], order[3], order[2], order[1]]

         


            '''
            for p in pairs:
                a,b = p["nt1"], p["nt2"]

                if {a,b} == {nt1, nt2}:
                    if p.get("lw") == "cHW":
                        order = [order[0], order[3], order[2], order[1]]
                    break
             '''

        


        
    
            nt1, nt2, nt3, nt4 = order
    
            self.nucl_quad[nt1] = Nucleotide(nucl[nt1], used_nucl, tetr_no, onz, 0)
            self.nucl_quad[nt2] = Nucleotide(nucl[nt2], used_nucl, tetr_no, onz, 1)
            self.nucl_quad[nt3] = Nucleotide(nucl[nt3], used_nucl, tetr_no, onz, 2)
            self.nucl_quad[nt4] = Nucleotide(nucl[nt4], used_nucl, tetr_no, onz, 3)
    
            self.tetrads.append(order)
            tetr_no += 1

    def GetChainFirstLast(self):
        chains = {}
        for _, nucl in self.nucl_quad.items():
            if nucl.chain not in chains:
                chains[nucl.chain] = {"first": "",\
                                      "last": "", \
                                      "val": sys.maxsize }
        # Analyze 
        curr_name = ""
        curr_min = sys.maxsize
        for name, nucl in self.nucl_quad.items():
            if chains[nucl.chain]["val"] > nucl.index:
                chains[nucl.chain]["val"] = nucl.index
                chains[nucl.chain]["first"] = name

        # Full connected_from variables.
        for chain, data in chains.items():
            curr_name = data["first"]
            next_conn = self.nucl_quad[data["first"]].connected_to
            while next_conn != "":
                self.nucl_quad[next_conn].connected_from = curr_name
                curr_name = next_conn
                next_conn = self.nucl_quad[next_conn].connected_to
            data["last"] = curr_name

        return chains

    def __init__(self, structure, quadruplex_id, tetrad_id = -1):
        self.nucl_quad = {}
        self.tetrads = []
        self.PrepareNucleotides(structure, quadruplex_id, tetrad_id)
        self.chains = self.GetChainFirstLast()
        if tetrad_id >= 0:
            single_tracts = []
            single_tracts.append(structure.tracts[quadruplex_id][tetrad_id])
            self.tracts = single_tracts

            # Loop
            single_loops = []
            single_loops.append(structure.loops[quadruplex_id][tetrad_id])
            self.loops = single_loops

        
        else:
            self.tracts = structure.tracts[quadruplex_id]
            self.loops = structure.loops[quadruplex_id]


        if os.getenv("draw_bulgeloop"):
            self.AddLoopBulge(structure)


    def GetNucleotidesPositions(self):
        lst = [-1] * len(self.nucl_quad)
        for name, nucl in self.nucl_quad.items():
            if nucl.connected_to != "":
              conn_to = self.nucl_quad[nucl.connected_to]
              lst[nucl.tetrade_no * 4 + nucl.position] = conn_to.tetrade_no * 4 + conn_to.position
            else:
              lst[nucl.tetrade_no * 4 + nucl.position] = -1
        return lst

    def GetSameRotations(self):
        lst = [-1] * len(self.tetrads)
        for idx, tetrad in enumerate(self.tetrads):
            tract_of = {}
            for group, tetrad_tracts in enumerate(self.tracts):
                for tract_index, names in enumerate(tetrad_tracts):
                    for nucl_name in names:
                        tract_of[nucl_name] = (group, tract_index)

            groups = set()
            for name in tetrad:
                if name in tract_of:
                    groups.add(tract_of[name][0])

            if len(groups) != 1:
                continue

            group = groups.pop()
            indices = [tract_of[name][1] for name in tetrad if name in tract_of]

            # Tetrads with the same tract direction can be aligned with each
            # other. A polarity-inverted block (e.g. in two-block quadruplexes)
            # has the opposite direction and must be aligned separately.
            # Otherwise the optimizer would have to reflect the tetrad, which
            # would flip its cWH/cHW pairing.
            direction = 0
            if len(indices) == 4:
                diffs = [(indices[(i + 1) % 4] - indices[i]) % 4 for i in range(4)]
                if all(diff == 1 for diff in diffs):
                    direction = 0
                elif all(diff == 3 for diff in diffs):
                    direction = 1

            lst[idx] = group * 2 + direction
        return lst

    def GetAlignments(self):
        lst = [-1] * len(self.nucl_quad)
        for name, nucl in self.nucl_quad.items():
            for group_1, tetrad_tracts in enumerate(self.tracts):
                for group_2, names in enumerate(tetrad_tracts):
                    for nucl_name in names:
                        if name == nucl_name:
                            lst[nucl.tetrade_no * 4 + nucl.position] = group_1 * 4 + group_2
        return lst

    # Use C++ code to rotate tetrads for more readable output.
    def Optimize(self, optimizer = "./svg_optimizer"):
        import optimizer
        optimized = optimizer.solve(self.GetNucleotidesPositions(),
                                    self.GetSameRotations(),
                                    self.GetAlignments())

        # Update position for nucleotide.
        for _, nucl in self.nucl_quad.items():
            pos = nucl.position + nucl.tetrade_no * 4
            for x in range(4):
                if int(optimized[nucl.tetrade_no * 4 + x]) == nucl.position:
                    nucl.position = x
                    break

     
        

        # Update positions in tetrades. For Tetrade border
        level = 0
        for tetrad in self.tetrads:
            new_positions = []
            new_positions.append(tetrad[int(optimized[level * 4 + 0])])
            new_positions.append(tetrad[int(optimized[level * 4 + 1])])
            new_positions.append(tetrad[int(optimized[level * 4 + 2])])
            new_positions.append(tetrad[int(optimized[level * 4 + 3])])

            #print("{0} -> {1}".format(tetrad, new_positions))
            for i, val in enumerate(new_positions):
                tetrad[i] = new_positions[i]

            level = level + 1


        if os.getenv("draw_5prime"):
           self.ForceFivePrime()

           for nucl in self.nucl_quad.values():
               nucl.connection_type = ConnType.UNKNOWN
               nucl.flow_in = ConnFlow.UNKNOWN
               nucl.flow_out = ConnFlow.UNKNOWN
               nucl.connected_from = ""

            
               nucl.priority_conn = -1
               nucl.priority_edge = -1
               nucl.priority_nucl = -1
                    
           self.chains = self.GetChainFirstLast()

           self.DetermineConnectionTypes()
           self.CalculateFlow(self.chains[next(iter(self.chains))])    

    
    

    def DetermineConnectionTypes(self):
        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "":
                conn = self.nucl_quad[nucl.connected_to]
                # Tetrades are counted from bottom to top.
                # > 0 - connection from top to bottom
                # < 0 - connection from bottom to top
                # = 0 - COnnection on the same level
                level_difference = conn.tetrade_no - nucl.tetrade_no
                if conn.position == nucl.position and abs(level_difference) == 1:
                    nucl.connection_type = ConnType.SIMPLE
                elif level_difference == 0:
                    nucl.connection_type = ConnType.SAME_LEVEL
                elif (nucl.position == 2 and conn.position == 2) or \
                     (nucl.position == 3 and conn.position == 3):
                    nucl.connection_type = ConnType.RIGHT
                elif (nucl.position == 2 and conn.position == 3) or \
                     (nucl.position == 3 and conn.position == 2):
                    nucl.connection_type = ConnType.RIGHT_CROSS
                elif (nucl.position == 0 and conn.position == 0) or \
                     (nucl.position == 1 and conn.position == 1):
                    nucl.connection_type = ConnType.LEFT
                elif (nucl.position == 0 and conn.position == 1) or \
                     (nucl.position == 1 and conn.position == 0):
                    nucl.connection_type = ConnType.LEFT_CROSS
                elif (nucl.position == 0 and conn.position == 3) or \
                     (nucl.position == 3 and conn.position == 0) or \
                     (nucl.position == 1 and conn.position == 2) or \
                     (nucl.position == 2 and conn.position == 1):
                    nucl.connection_type = ConnType.FRONT_BACK_CROSS
                elif (nucl.position == 1 and conn.position == 3) or \
                     (nucl.position == 3 and conn.position == 1) or \
                     (nucl.position == 0 and conn.position == 2) or \
                     (nucl.position == 2 and conn.position == 0):
                    nucl.connection_type = ConnType.FRONT_TO_BACK


    def PrintFlow(self, chain):
        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            #print(nucl.full_name, nucl.flow_out, nucl.flow_in)
            nucl = conn

    # TODO Optimize. Do I rly need so many iterations?
    # I HAVE TO FIX IT!! it "works"
    def CalculateFlow(self, chain):
        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "":
                conn = self.nucl_quad[nucl.connected_to]
                # Tetrades are counted from bottom to top.
                # > 0 - connection from top to bottom
                # < 0 - connection from bottom to top
                # = 0 - COnnection on the same level
                level_difference = nucl.tetrade_no - conn.tetrade_no
                if conn.position == nucl.position and abs(level_difference) == 1:
                    if level_difference == -1:
                        # Connection going bot to top
                        nucl.flow_out = ConnFlow.UP
                        nucl.flow_in = ConnFlow.UP
                        conn.flow_out = ConnFlow.UP
                    else:
                        # Connection top to bot
                        nucl.flow_out = ConnFlow.DOWN
                        nucl.flow_in = ConnFlow.DOWN
                        conn.flow_out = ConnFlow.DOWN


        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                conn = self.nucl_quad[nucl.connected_to]
                if conn.flow_out == ConnFlow.UP:
                    nucl.flow_in = ConnFlow.DOWN
                elif conn.flow_out == ConnFlow.DOWN:
                    nucl.flow_in = ConnFlow.UP

        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]

            if nucl.flow_in == ConnFlow.UNKNOWN:
                level_difference = nucl.tetrade_no - conn.tetrade_no
                if level_difference < 0:
                    # Connection going bot to top
                    nucl.flow_in = ConnFlow.DOWN
                    conn.flow_out = ConnFlow.UP
                elif level_difference > 0:
                    # Connection top to bot
                    nucl.flow_in = ConnFlow.UP
                    conn.flow_out = ConnFlow.DOWN
            nucl = conn

        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                conn = self.nucl_quad[nucl.connected_to]
                if conn.flow_out == ConnFlow.UP:
                    nucl.flow_in = ConnFlow.DOWN
                elif conn.flow_out == ConnFlow.DOWN:
                    nucl.flow_in = ConnFlow.UP

        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            if nucl.flow_in == ConnFlow.UNKNOWN:
                nucl.flow_in = nucl.flow_out
                if nucl.flow_out == ConnFlow.UP:
                    conn.flow_out = ConnFlow.DOWN
                elif nucl.flow_out == ConnFlow.DOWN:
                    conn.flow_out = ConnFlow.UP
            nucl = conn

        nucl = self.nucl_quad[chain["last"]]
        while nucl.connected_from != "":
            conn = self.nucl_quad[nucl.connected_from]
            if conn.flow_out == ConnFlow.UNKNOWN:
                level_difference = nucl.tetrade_no - conn.tetrade_no
                if level_difference < 0:
                    # Connection going bot to top
                    conn.flow_out = ConnFlow.DOWN
                elif level_difference > 0:
                    # Connection top to bot
                    conn.flow_out = ConnFlow.UP
                else:
                    conn.flow_out = conn.flow_in
            nucl = conn

        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                conn = self.nucl_quad[nucl.connected_to]
                if conn.flow_out == ConnFlow.UP:
                    nucl.flow_in = ConnFlow.DOWN
                elif conn.flow_out == ConnFlow.DOWN:
                    nucl.flow_in = ConnFlow.UP

        # Determine starting flow.
        nucl = self.nucl_quad[chain["first"]]
        if nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            level_difference = nucl.tetrade_no - conn.tetrade_no
            if level_difference < 0:
                # Connection going bot to top
                nucl.flow_out = ConnFlow.UP
            elif level_difference > 0:
                # Connection top to bot
                nucl.flow_out = ConnFlow.DOWN
            else:
                nucl.flow_out = nucl.flow_in

            for _, nucl in self.nucl_quad.items():
                # Check if it is not the last nucleotide.
                if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                    conn = self.nucl_quad[nucl.connected_to]
                    if conn.flow_out == ConnFlow.UP:
                        nucl.flow_in = ConnFlow.DOWN
                    elif conn.flow_out == ConnFlow.DOWN:
                        nucl.flow_in = ConnFlow.UP

        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            if nucl.flow_in == ConnFlow.UNKNOWN:
                if nucl.flow_out == ConnFlow.UNKNOWN:
                    if nucl.connected_from == "":
                        nucl.flow_out = ConnFlow.UP
                    else:
                        conn_from = self.nucl_quad[nucl.connected_from]
                        if conn_from.flow_in == ConnFlow.DOWN:
                            nucl.flow_out = ConnFlow.UP
                        else:
                            nucl.flow_out = ConnFlow.DOWN
                nucl.flow_in = nucl.flow_out
                if nucl.flow_out == ConnFlow.UP:
                    conn.flow_out = ConnFlow.DOWN
                elif nucl.flow_out == ConnFlow.DOWN:
                    conn.flow_out = ConnFlow.UP
            nucl = conn

        #print("\n\n")
        #self.PrintFlow(chain)

class Structure:
    def __init__(self):
        self.nucleotides = {}
        # Index represents corepsponding quadruplex ID data.
        self.tetrads = []
        self.tetrads_order = []
        self.single_tetrads = []
        self.tracts = []
        self.loops = []



    """
    data - json dict of tetrade parameters.
    {
      "id":"A.DC2-A.DG10-B.DC9-B.DG3",
      "nt1":"A.DC2",
      "nt2":"A.DG10",
      "nt3":"B.DC9",
      "nt4":"B.DG3",
      "onz":"N+",
      "gbaClassification":"VIIIa",
      "planarityDeviation":0.1270031987786134,
      "ionsChannel":[
      ],
      "ionsOutside":[
      ]
    },
    """
    #def addTetrade(self, name, data, tetrade_no):
    #    self.tetrads[tetrade_no][name] = data

    """
    data - json dict of nucleotide parameters.
    {
      "index":22,
      "model":1,
      "chain":"B",
      "number":11,
      "icode":" ",
      "molecule":"DNA",
      "fullName":"B.DG11",
      "shortName":"G",
      "chi":-89.97,
      "glycosidicBond":"syn"
    }
    """
    def addNucleotide(self, name, data):
        self.nucleotides[name] = data

    def fromFile(self, path):
        with open(path) as file:
            data = json.load(file)
        return self.fromJsonDict(data)

    def fromString(self, json_string):
        data = json.loads(json_string)
        return self.fromJsonDict(data)

    def fromJsonDict(self, json_dict):
        for data in json_dict["nucleotides"]:
            self.addNucleotide(data["fullName"], data)
            
        self.basePairs = json_dict.get("basePairs", [])


        # TODO What to do with other helices/quadruplexes?
        for index, helice in enumerate(json_dict["helices"]): #["quadruplexes"][0]["tetrads"].items():
            # Unordered tetrads. Order them from "tetrad_pairs"
            single_tetrads_local = []
            tetrad_unordered = {}
            tracts_all = []
            loops_all = []
            for _, quadruplex in enumerate(helice["quadruplexes"]):
                tetrad_unordered_local = {}
                for data in quadruplex["tetrads"]:
                    tetrad_unordered[data["id"]] = data
                    tetrad_unordered_local[data["id"]] = data
                single_tetrads_local.append(tetrad_unordered_local)
                if "tracts" in quadruplex:
                    tracts_all.append(quadruplex["tracts"])
                else:
                    tracts_all.append(list())
                    
                if "loops" in quadruplex:
                    loops_all.append(quadruplex["loops"])
                else:
                    loops_all.append(list())

            # Order tetrads according to "tetrad_pairs" data.
            tetrad_pairs = helice["tetradPairs"]
            tetrad_ordered = []
            if tetrad_pairs != None and len(tetrad_pairs) > 0:
                # Add first pair.
                pair = tetrad_pairs.pop(0)
                tetrad_ordered.append(pair["tetrad1"])
                tetrad_ordered.append(pair["tetrad2"])
                to_remove = 0
                while len(tetrad_pairs) > 0:
                    found = False
                    for index_pair, pair in enumerate(tetrad_pairs):
                        for index, tetrad in enumerate(tetrad_ordered):
                            if tetrad == pair["tetrad1"]:
                                tetrad_ordered.insert(index + 1, pair["tetrad2"])
                                found = True
                                break
                            elif tetrad == pair["tetrad2"]:
                                tetrad_ordered.insert(index, pair["tetrad1"])
                                found = True
                                break
                        if found == True:
                            to_remove = index_pair
                            break
                    tetrad_pairs.pop(to_remove)
                    if len(tetrad_pairs) == 0:
                        break

            tetrad_ordered.reverse()

            
            
            
            # Do not add single tetrads as quadruplexes.
            if len(tetrad_ordered) > 1:
                self.tetrads_order.append(tetrad_ordered)
                self.tetrads.append(tetrad_unordered)

            self.tracts.append(tracts_all)
            self.loops.append(loops_all)
            self.single_tetrads.append(single_tetrads_local)
        
        
        
        
        # TODO check if input data was valid.
        # Are there any nucleotides and tetrads?
        return self
