# Copyright (C) <2022>  <The Ohio State University>       

# This program is free software: you can redistribute it and/or modify                              
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or    
# (at your option) any later version.                                                                                       
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
# GNU General Public License for more details.                                                                             

# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

def fasta_to_dict(fasta_in):
    """
    """
    seqid_dict = {}
    with open(fasta_in, "r") as f:
        subj, seq = "", ""

        for row in f:
            if row.startswith(">"):
                if len(subj) != 0:
                    seqid_dict[subj] = seq
                subj = row[1:].replace("\n","").split()[0]
                seq = ""
            else:
                seq = seq + row.replace("\n","")
        
        if len(subj) != 0:
            seqid_dict[subj] = seq
    f.close()

    return seqid_dict

