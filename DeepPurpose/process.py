import pickle
from Bio import SeqIO

import pandas as pd
import numpy as np

# 6 unique drugs in BindingDB with KD
chem_bdb = {'Nc1nc(N)c2cc(NCc3ccc(cc3)[N+]([O-])=O)ccc2n1',
            'Nc1nc(N)c2cc(ccc2n1)[N+]([O-])=O',
            'Nc1nc(N)c2cc(NCc3ccc(O)cc3)ccc2n1',
            'COc1ccc(CNc2ccc3nc(N)nc(N)c3c2)cc1',
            'Nc1nc(N)c2cc(NCCn3c(nc4cc(Cl)c(Cl)cc34)C(F)(F)F)ccc2n1',
            'Nc1ccc2nc(N)nc(N)c2c1'}

# 2 unique protein in BindingDB with KD
prot_bdb = {'MTAPTVPVALVTGAAKRLGRSIAEGLHAEGYAVCLHYHRSAAEANALSATLNARRPNSAITVQADLSNVATAPVSGADGSAPVTLFTRCAELVAACYTHWGRCDVLVNNASSFYPTPLLRNDEDGHEPCVGDREAMETATADLFGSNAIAPYFLIKAFAHRFAGTPAKHRGTNYSIINMVDAMTNQPLLGYTIYTMAKGALEGLTRSAALELAPLQIRVNGVGPGLSVLVDDMPPAVWEGHRSKVPLYQRDSSAAEVSDVVIFLCSSKAKYITGTCVKVDGGYSLTRA',
            'MSRAAARFKIPMPETKADFAFPSLRAFSIVVALDMQHGIGDGESIPWRVPEDMTFFKNQTTLLRNKKPPTEKKRNAVVMGRKTWESVPVKFRPLKGRLNIVLSSKATVEELLAPLPEGQRAAAAQDVVVVNGGLAEALRLLARPLYCSSIETAYCVGGAQVYADAMLSPCIEKLQEVYLTRIYATAPACTRFFPFPPENAATAWDLASSQGRRKSEAEGLEFEICKYVPRNHEERQYLELIDRIMKTGIVKEDRTGVGTISLFGAQMRFSLRDNRLPLLTTKRVFWRGVCEELLWFLRGETSAQLLADKDIHIWDGNGSREFLDSRGLTENKEMDLGPVYGFQWRHFGADYKGFEANYDGEGVDQIKLIVETIKTNPNDRRLLVTAWNPCALQKMALPPCHLLAQFYVNTDTSELSCMLYQRSCDMGLGVPFNIASYALLTILIAKATGLRPGELVHTLGDAHVYRNHVDALKAQLERVPHAFPTLIFKEERQYLEDYELTDMEVIDYVPHPAIKMEMAV'}

# proteins
path = "data/Targets/l.major.fasta"
records = list(SeqIO.parse(path, "fasta"))
prot_cha = {i.seq._data for i in records}               # 8,495 proteins

path = "data/Targets/preferredTargets.unique.fasta"
records_pref = list(SeqIO.parse(path, "fasta"))
prot_cha_pref = {i.seq._data for i in records_pref}     # 34,594 proteins

path = "data/Targets/all_targets.fasta"
records_all = list(SeqIO.parse(path, "fasta"))
prot_cha_all = {i.seq._data for i in records_all}       # 79,982 proteins

# chemicals
ddd = pd.read_csv("data/Molecules/drugBank_leishmania.smiles")
d = set(ddd['smiles'])

ddd = pd.read_csv("data/Molecules/drugCentral.csv", sep=',', header=0, usecols=[1, 4])
dp = set(ddd['SMILES'])

ddd = pd.read_csv("data/Molecules/endogenous.csv", sep=',', header=0)
dpp = set(ddd['smiles'])

ddd = pd.read_csv("data/Molecules/in-trials.csv", sep=',', header=0)
dppp = set(ddd['smiles'])

ddd = pd.read_csv("data/Molecules/world.csv", sep=',', header=0)
dpppp = set(ddd['smiles'])

chem_all = d | dp | dpp | dppp | dpppp                  # 94,053 chemiclas

# proteins
overlap_prot = prot_cha & prot_bdb              # 2 overlapping protein with l.major
overlap_prot_pref = prot_cha_pref & prot_bdb    # 2 overlapping protein with l.prefferredTargets
overlap_prot_all = prot_cha_all & prot_bdb      # 2 overlapping protein with l.all

# chemicals
overlap_chem = chem_all & chem_bdb              # No overlapping chemiclas

from rdkit import Chem
from rdkit.DataStructs import FingerprintSimilarity
m1 = Chem.MolFromSmiles('Nc1ccc2nc(N)nc(N)c2c1')
m2 = Chem.MolFromSmiles('COc1ccc(CNc2ccc3nc(N)nc(N)c3c2)cc1')
FingerprintSimilarity(Chem.RDKFingerprint(m1), Chem.RDKFingerprint(m2))

n_bdb = len(chem_bdb)
chems = list(chem_bdb).extend(list(chem_all))


import requests
pd.DataFrame(data=[list(chem_all), ['CHEM_{}'.format(i) for i in range(1, len(chem_all+1))]])


# ##########################################
# BLAST AAs
from Bio.PDB import parse_pdb_header
import json
import requests


def blast_aa_seq(aa_sequence, i):
    page = """https://search.rcsb.org/rcsbsearch/v1/query?json=
              {"query": {"type": "terminal", "service": "sequence",
                         "parameters": {"evalue_cutoff": 1, 
                                       "identity_cutoff": 1, 
                                       "target": "pdb_protein_sequence",
                                       "value": \"""" + str(aa_sequence) + """\"}},
               "request_options": {"scoring_strategy": "sequence"},
               "return_type": "polymer_entity"}"""

    req = requests.get(page)
    if req.status_code == 200:
        return req.text
    else:
        print("no response - {}".format(i))
        return None


for i in range(len(records_all)):
    seq = records_all[i].seq._data
    response = blast_aa_seq(seq, i)
    pdb_reso = []

    try:
        results = json.loads(response)['result_set']
        for r in results:
            pdb_id = None
            pdb_full = None
            # match
            info = r['services'][0]['nodes'][0]['match_context'][0]
            if info['mismatches'] == 0 and info['gaps_opened'] == 0 and info['query_length'] == info['subject_length']:
                pdb_id = r['identifier'].split('_')[0]
                pdb_full = r['identifier']

            # if match, download pdb file
            if pdb_id and pdb_full:
                page = 'http://files.rcsb.org/view/{}.pdb'.format(pdb_id)
                req = requests.get(page)
                if req.status_code == 200:
                    response = req.text
                    outfile = 'tmp.pdb'
                    if outfile:
                        with open(outfile, 'w') as f:
                            f.write(response)
                        # parse to get the resolution
                        structure = parse_pdb_header(outfile)
                        pdb_reso.append((pdb_full, structure['resolution']))

        # append to dataset file
        if pdb_reso:
            # find the pdb with best resolution
            tmp_dict = {r: p for p, r in pdb_reso}
            best_pdb_id = tmp_dict[max(tmp_dict.keys())]
            # write to file
            with open('./all_targets.csv', 'a') as f:
                f.write("{}, {}, {}, {}\n".format(best_pdb_id, pdb_reso, records_all[i].description, seq))
                print("{} - {}".format(i, pdb_reso))
    except:
        pass

# #################################################################
# #################################################################
# 2021/02/02
# evaluate results' rank correlation between DeepPurpose and MONN
# #################################################################
# #################################################################

import pandas as pd
protein_deep = pd.read_csv('data/pdb_leish_all.tsv', sep='\t', usecols=[0, 4])

protein_table = pd.DataFrame(
    [['Q6TUJ5', '4iu0', 'MEHVQQYKFYKEKKMSIVLAPFSGGQPHSGVELGPDYLLKQGLQQDMEKLGWDTRLERVFDGKVVEARKASDNGDRIGRVKRPRLTAECTEKIYKCVRRVAEQGRFPLTIGGDHSIALGTVAGVLSVHPDAGVIWVDAHADINTMSGTVSGNLHGCPLSILLGLDRENIPECFSWVPQVLKPNKIAYIGLRAVDDEEKKILHDLNIAAFSMHHVDRYGIDKVVSMAIEAVSPKGTEPVMVSYDVDTIDPLYVPATGTPVRGGLSFREALFLCERIAECGRLVALDVVECNPLLAATESHVNDTISVGCAIARCMMGETLLYTPHTSSKL'],
     ['P48499', '1amk', 'MSAKPQPIAAANWKCNGTTASIEKLVQVFNEHTISHDVQCVVAPTFVHIPLVQAKLRNPKYVISAENAIAKSGAFTGEVSMPILKDIGVHWVILGHSERRTYYGETDEIVAQKVSEACKQGFMVIACIGETLQQREANQTAKVVLSQTSAIAAKLTKDAWNQVVLAYEPVWAIGTGKVATPEQAQEVHLLLRKWVSENIGTDVAAKLRILYGGSVNAANAATLYAKPDINGFLVGGASLKPEFRDIIDATR'],
     ['O97193', '5ofu', 'MDVRRTPTPTTLTQYIIKSQPPHSRGDFTLLMMAIQTSVKVIEKNIRRAGMKGMLGYIAGQSANATGDHQAKLDVISNIAFKAYLLSSTSVCVLGSEEEEQMIIAESGRRGDYLIFFDPLDGSSNIDANVSVGSIWGVWRLPKDTTINSVEDANAVIRMLKGTDMVSAGYAVYGSATNLVLTSGHGVDGFTLDPNIGEFILTHPHISIPKKRSIYSVNEGNYGKWEPWFKEYIDYLKMNKTTRYSARYIGSMVGDIHRTLLYGGIFCYPKDANQVEGKLRLLYEAAPMAMIVEQAGGKAVGSNGRILEQSITRLHQRTPVYFGSRQEVDLCMAFRDRNVKTEALAPTSSKL'],
     ['O15826', '2yay', 'MKRARSANIPGAILHSLAELQDGLNAMIDPSWRAVRSLDNWALAITMESTELLDSYPWKWWKNLNATPDLANVRIELVDIFHFSLSGAMQMRSTPDDEIPAASLKPLKEVMTTFLPAKECTSDPYGFVFFPLTDTQNAIASFRNIIQLANAYRFDVIIECIIYAAEDLGFNLVAYYIAKHTLNCIRQLSGYKDGSYVKVNNGVEDNSLLHNCIKDVSLDEVLDADKYVQAWNSIMANVYEAFQIKESDRKDAERWFALAKENRLAIKA'],
     ['Q27686', '3pp7', 'MSQLAHNLTLSIFDPVANYRAARIICTIGPSTQSVEALKGLIQSGMSVARMNFSHGSHEYHQTTINNVRQAAAELGVNIAIALDTKGPEIRTGQFVGGDAVMERGATCYVTTDPAFADKGTKDKFYIDYQNLSKVVRPGNYIYIDDGILILQVQSHEDEQTLECTVTNSHTISDRRGVNLPGCDVDLPAVSAKDRVDLQFGVEQGVDMIFASFIRSAEQVGDVRKALGPKGRDIMIICKIENHQGVQNIDSIIEESDGIMVARGDLGVEIPAEKVVVAQKILISKCNVAGKPVICATQMLESMTYNPRPTRAEVSDVANAVFNGADCVMLSGETAKGKYPNEVVQYMARICLEAQSALNEYVFFNSIKKLQHIPMSADEAVCSSAVNSVYETKAKAMVVLSNTGRSARLVAKYRPNCPIVCVTTRLQTCRQLNITQGVESVFFDADKLGHDEGKEHRVAAGVEFAKSKGYVQTGDYCVVIHADHKVKGYANQTRILLVE'],
     ['Q4QHJ8', '3uib', 'MQAKGEAAMRDLIAELHAMQSPYTVQRFISSGSYGAVCAGVDSEGIPVAIKRVFNTVSDGRTVNILSDSFLCKRVLREIRLLNHFHHPNILGLRDIFVHFEEPAMHKLYLVTELMRTDLAQVIHDQRIVISPQHIQYFMYHILLGLHVLHEAGVVHRDLHPGNILLADNNDITICDFNLAREDTADANKTHYVTHRWYRAPELVMQFKGFTKLVDMWSAGCVMAEMFNRKALFRGSTFYNQLNKIVEVVGTPKIEDVVMFSSPSARDYLRNSLSNVPARAWTAVVPTADPVALDLIAKMLEFNPQRRISTEQALRHPYFESLFDPLDLTEGLSERFHFDESVTDVYDMHKIFTAEVERFNDLRERREEVARERAVAAQQQGEQVLGTDHMPRTHSLMELAGSAPAPS'],
     ['Q4QBL1', '4k10', 'MAHMERFQKVYEEVQEFLLGDAEKRFEMDVHRKGYLKSMMDTTCLGGKYNRGLCVVDVAEAMAKDTQMDAAAMERVLHDACVCGWMIEMLQAHFLVEDDIMDHSKTRRGKPCWYLHPGVTAQVAINDGLILLAWATQMALHYFADRPFLAEVLRVFHDVDLTTTIGQLYDVTSMVDSAKLDAKVAHANTTDYVEYTPFNHRRIVVYKTAYYTYWLPLVMGLLVSGTLEKVDKKATHKVAMVMGEYFQVQDDVMDCFTPPEKLGKIGTDIEDAKCSWLAVTFLTTAPAEKVAEFKANYGSTDPAAVAVIKQLYTEQNLLARFEEYEKAVVAEVEQLIAALEAQNAAFAASVKVLWSKTYKRQK']],
    columns=['uniprot_id', 'pdb_id', 'aas']
)

dp_set = set([i[1:] for i in protein_deep['AAs'].tolist()])
for i, j, k in protein_table.to_numpy():
    if k in dp_set:
        print(i, k)
# --> there are only two overlapping proteins: P48499 and Q27686

# For the in-trials datasets
monn_result = pd.read_csv('MONN/preds_in-trial_ordered.csv')
# chem_input_table = pd.DataFrame({'smiles': ddd['smiles'], 'index': ddd['zinc_id']})

chem_input_table = pd.DataFrame({'smiles': ddd['SMILES'], 'index': ['id-'+str(i) for i in ddd['ID']]})
# chem_input_table = pd.DataFrame({'smiles': ddd['smiles'], 'index': ['id-'+str(i) for i in range(ddd.shape[0])]})
chem_input_table.to_csv('MONN/data/drugcentral.csv', index=False)

# merge results on proteins for each drug dataset
drugset = ['in-trial', 'drugbank', 'drugcentral', 'endogenous', 'world']
protset = ['Q6TUJ5', 'O97193', 'O15826', 'Q4QHJ8', 'Q4QBL1', 'P48499', 'Q27686']

i = 0
j = 0

for i in range(len(drugset)):
    out = []
    for j in range(7):
        with open(
                'MONN/{}/save_folder_{}/results_aggregation/output_list.pkl'.format(
                        drugset[i], protset[j]), 'rb') as f:
            out.append(pickle.load(f))
    out = pd.DataFrame(np.array(out).reshape((-1, 3)),
                       columns=['drug', 'prot', 'score'])
    out.sort_values(by=['score'], inplace=True)
    out['score'] = out['score'].astype(np.float)
    out.to_csv('MONN/DeepPurpose_{}_rank.csv'.format(drugset[i]), index=False)
    a = out['score'][out['drug'] == 'ZINC000022059930'][out['prot'] == 'Q6TUJ5']

    # compute correlation
    monn = pd.read_csv('MONN/preds_{}_ordered.csv'.format(drugset[i]),
                       usecols=[0, 1, 3])
    f = open('MONN/{}/merged_rank.csv'.format(drugset[i]), 'w')
    f.write("zinc_id,pid,Kd,DeepPurpose_score\n")
    for z, p, s in monn.to_numpy():
        tmp = out['score'][out['drug'] == z][out['prot'] == p].to_numpy()
        tmp = tmp[0] if tmp.tolist() else np.NAN
        f.write("{},{},{},{}\n".format(z, p, s, tmp))
    f.close()

from scipy.stats import spearmanr
result_all = pd.read_csv('MONN/in-trial/in-trails_merged_rank.csv', usecols=[2, 3])
result_all = result_all.dropna()
correlation, pvalue = spearmanr(result_all['Kd'].tolist(), result_all['DeepPurpose_score'].tolist())

