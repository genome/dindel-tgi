import sys, os, Fasta
class Variant:
    def __init__(self, varString = ''):
        
        l = len(varString)

        if varString[0] == '-' and l>1:
            self.type = "del"
            self.seq = varString[1:]
            self.length = len(varString)-1
        elif varString[0] == '+' and l>1:
            self.type = "ins"
            self.seq = varString[1:]
            self.length = len(varString)-1
        elif len(varString)==4 and varString[1:3]=='=>':
            self.type = "snp"
            self.seq = varString[3]
            self.length = 1
        elif varString[0] == '*' or varString.find("REF") != -1 or varString.find("ref") != -1:
            self.type = "ref"
            self.length = 0
            self.seq = ''
        else:
            raise NameError("Unrecognized variant: "+varString)
        
        # offset is used when the variant is initialized using an alt/ref pair
        self.offset = 0

        self.str = varString

class Variant4:
    def __init__(self, ref = '', alt = ''):
        # create Dindel style variant from ref/alt combination
        
        dlen = len(ref) - len(alt)
        
        if dlen == 0:
            nm = 0
            altnuc = ''
            refnuc = ''
            for idx in range(len(ref)):
                a = ref[idx]
                b = alt[idx]
                if a!=b:
                    nm += 1
                    self.offset = idx
                    refnuc = a
                    altnuc = b


            if nm == 0:
                self.type = "ref"
                self.length = 0
                self.seq = ''
                self.str = 'REF'
            elif nm == 1:
                self.type = "snp"
                self.length = 1
                self.seq = altnuc
                self.str = "%s=>%s" % (refnuc, altnuc)
            else:
                raise NameError("MultiSNP")
        else:
            # AT
            # ACCCT

            if dlen<0:
                self.type = "ins"
                _alt = alt[:]
                _ref = ref[:]
                self.str = '+'
            else:
                self.type = "del"
                _alt = ref[:]
                _ref = alt[:]
                self.str = '-'

            numrb = len(_ref) # number of reference bases in alt that need to be identified 
            left_match = 0
            right_match = 0
            for x in range(0, len(_ref)+1):
                if _ref[:x] == _alt[:x]:
                    left_match = x
            for x in range(0, len(_ref)+1):
                if _ref[-x:] == _alt[-x:]:
                    right_match = x

            if left_match == 0 or left_match + right_match < numrb:
                raise NameError("Don't think this is a proper VCF4 insertion")
            
            # try to shift it to the left as much as possible
            left_end = 1

            if numrb-left_end>right_match:
                left_end = left_match

            right_start = numrb-left_end
            if right_start == 0:
                right_start = -len(_alt)
            self.seq = _alt[left_end:-right_start]    
            self.offset = left_end
            self.str += self.seq
            self.length = len(self.seq)
           
    
            #print self.type, self.seq, self.offset, left_end

            





def isIndel(allele =''):
    if allele.find('/')!=-1:
        raise NameError('Is not allele, but genotype!')
    if allele[0]=='+' or allele[0]=='-':
        return True
    else:
        return False
