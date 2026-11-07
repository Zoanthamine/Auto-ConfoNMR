"""Exercise real generated stage inputs and unchanged noniodine behavior."""
import copy
import unittest
from iodine_compat import configure,render,require_policy


class IodineCompatibilityChecks(unittest.TestCase):
    def cfg(self):return {'max_conformers':6,'postprocessing':[{'allow_exploratory':True,'differences':[]}]}
    def template(self,route):return '%mem=60GB\n'+route+'\n\nname\n\n0 1\nC 0 0 0\nI 2.1 0 0\n\n'

    def test_ordinary_input_and_config_are_byte_identical(self):
        cfg=self.cfg();before=copy.deepcopy(cfg);text=self.template('# B3LYP/6-31G* opt')
        self.assertEqual(configure(cfg,['C','H']),before)
        self.assertEqual(render(text,['C','H'],cfg,'opt'),text)

    def test_unsupported_iodine_stops_before_calculation(self):
        with self.assertRaisesRegex(ValueError,'default 6-31G/6-311G'):require_policy(self.cfg(),['C','I'])

    def test_opt_and_resp_pair_iodine_basis_with_ecp(self):
        cfg=configure(self.cfg(),['C','H','I'])
        for stage in ['resp','opt','freq']:
            route='# B3LYP/6-31G* '+('pop=MK' if stage=='resp' else 'opt')
            text=render(self.template(route),['C','H','I'],cfg,stage)
            self.assertIn('B3LYP/GenECP',text)
            self.assertIn('C H 0\n6-31G*\n****',text)
            self.assertIn('I 0\nDef2SVP\n****\n\nI 0\nDef2SVP',text)
            self.assertIn('6D',text)
            if stage=='resp':
                self.assertIn('Pop=(MK,ReadRadii)',text)
                self.assertTrue(text.endswith('I 1.98\n\n'))

    def test_nmr_keeps_functional_light_basis_and_solvent_uses_all_electron_i(self):
        cfg=configure(self.cfg(),['C','H','I'])
        route='# mPW1PW91/6-311+G(d,p) NMR=GIAO SCRF=(IEFPCM,Solvent=Chloroform) G09Defaults'
        text=render(self.template(route),['C','H','I'],cfg,'nmr')
        self.assertIn('mPW1PW91/Gen NMR=GIAO SCRF=(IEFPCM,Solvent=Chloroform) G09Defaults',text)
        self.assertIn('C H 0\n6-311+G(d,p)',text)
        self.assertIn('I 0\nDGDZVP',text)
        self.assertNotIn('GenECP',text)
        self.assertFalse(cfg['postprocessing'][0]['allow_exploratory'])
        self.assertEqual(cfg['max_conformers'],6)

    def test_custom_basis_not_overwritten(self):
        cfg=configure(self.cfg(),['C','I'])
        with self.assertRaisesRegex(ValueError,'custom inputs require review'):
            render(self.template('# B3LYP/Gen'),['C','I'],cfg,'opt')


if __name__=='__main__':unittest.main()
