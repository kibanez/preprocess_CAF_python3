import json
import os
import numpy as np
import pandas as pd
import epydemiology as epy
from datetime import datetime


class CAFPreprocess:

    def __init__(self, hash_cfg):
        # Input parameters
        self.id_cases = hash_cfg.get('id_cases')
        if self.id_cases not in ("CGRSequenceID", "SampleID"):
            raise IOError('The ID used for the list of cases is missing - '
                          'Please specify `CGRSequenceID` or `SampleID`')
        f_cases = hash_cfg.get('list_cases')
        if not os.path.isfile(f_cases):
            raise IOError('The file containing the list of cases %s dose not exist'
                          % f_cases)
        else:
            self.l_cases = self.read_cohort(f_cases)

        self.id_controls = hash_cfg.get('id_controls')
        if self.id_controls not in ("CGRSequenceID", "SampleID"):
            raise IOError('The ID used for the list of controls is missing - '
                          'Please specify `CGRSequenceID` or `SampleID`')
        f_controls = hash_cfg.get('list_controls')
        if not os.path.isfile(f_controls):
            raise IOError('The file containing the list of controls %s dose not exist'
                          % f_controls)
        else:
            self.l_controls = self.read_cohort(f_controls)

        # matched cases and control lists
        self.l_cases_matched = []
        self.l_controls_matched = []

        f_metadata = hash_cfg.get('metadata_file')
        if not os.path.isfile(f_metadata):
            raise IOError('The QC metadata file %s does not exist' % f_metadata)
        else:
            self.metadata_cases = self.read_metadata(f_metadata,
                                                     self.l_cases,
                                                     self.id_cases)
            self.metadata_controls = self.read_metadata(f_metadata,
                                                        self.l_controls,
                                                        self.id_controls)

        verify_bam_id = hash_cfg.get('verifybamid')
        if verify_bam_id is None:
            raise IOError('An integer is expected for the contamination threshold')
        else:
            self.verify_bam_id = float(verify_bam_id)

        cc_ratio = hash_cfg.get('cc_ratio')
        if cc_ratio is None:
            raise IOError('A ratio for cases-control is required')
        else:
            self.cc_ratio = cc_ratio

        match_gender = hash_cfg.get('match_gender')
        if match_gender is None:
            raise IOError('A boolean value is expected for the gender match')
        else:
            self.match_gender = True if match_gender == 'True' else False

        match_age = hash_cfg.get('match_age')
        if match_age is None:
            raise IOError('A boolean value is expected for the age match')
        else:
            self.match_age = True if match_age == 'True' else False

        match_ancestry = hash_cfg.get('match_ancestry')
        if match_ancestry is None:
            raise IOError('A boolean value is expected for the ancestry match')
        else:
            self.match_ancestry = True if match_ancestry == 'True' else False

        # Output parameters
        output_path = hash_cfg.get('output_path')
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        self.output_path = output_path
        self.output_cases = os.path.join(self.output_path,
                                         os.path.basename(f_cases).split(".txt")[0] +
                                         "_" +
                                         datetime.today().strftime('%Y-%m-%d') +
                                         ".txt")
        self.output_controls = os.path.join(self.output_path,
                                            os.path.basename(f_controls).split(".txt")[0] +
                                            "_" +
                                            datetime.today().strftime('%Y-%m-%d') +
                                            ".txt")
        self.output_cases_matched = os.path.join(self.output_path,
                                         "list_cases_matched_" +
                                         datetime.today().strftime('%Y-%m-%d') + ".txt")
        self.output_controls_matched = os.path.join(self.output_path,
                                            "list_controls_matched_" +
                                            datetime.today().strftime('%Y-%m-%d') + ".txt")
        self.output_log_file = os.path.join(self.output_path,
                                            "preprocessing_CAF_output" +
                                            datetime.today().strftime('%Y-%m-%d') +
                                            ".json")
        # We init output dict first with the info we get from the config file
        self.output_log = hash_cfg

    def read_cohort(self, file_cohort):
        """
        Reads a file containing individuals within a cohort (cases/controls) and
        it returns a list of them
        :param file_cohort: path to the file containing the list of individuals
        :return: list of CGRSequenceID
        """
        cohort_data = pd.read_csv(file_cohort, sep='\t',
                                  error_bad_lines=False, header=None)
        return cohort_data.iloc[:, 0].tolist()

    def read_metadata(self, file_metadata, l_cohort, id_selection):
        """
        Reads a file with the QC metadata created by CGR internal bioinformatics pipeline
        :param file_metadata: path to the QC metadata file
        :param l_cohort: list of individuals within that cohort (cases or controls)
        :param id_selection: binary value: "CGRSequenceID" or "SampleID"
        :return: dict with all metadata columns
        """
        qc_metadata = pd.read_csv(file_metadata, sep=',', error_bad_lines=False)
        return qc_metadata.loc[qc_metadata[id_selection].isin(l_cohort)]

    def filter_out_contamination(self):
        """
        It filters out and removes contaminated individuals from both cases and control lists
        :return: nothing. updates self.l_controls and self.l_cases (and the log dic)
        """
        cases_cont = \
            self.metadata_cases.loc[self.metadata_cases['verifybamid_freemix'] >
                                    self.verify_bam_id]['CGRSequenceID'].tolist()
        controls_cont = \
            self.metadata_controls.loc[self.metadata_controls[
                                                              'verifybamid_freemix'] >
                                    self.verify_bam_id]['CGRSequenceID'].tolist()
        self.output_log['cases_contamination'] = cases_cont
        self.output_log['controls_contamination'] = controls_cont
        # Update the list of cases and controls
        self.l_cases = list(set(self.l_cases) - set(cases_cont))
        self.l_controls = list(set(self.l_controls) - set(controls_cont))

    def write_log(self):
        """
        Function that writes in self.output_log_file the summary of the whole process, contained
        in self.output_log dict. It also updates the lists of cases and controls from what we
        have filtered out (depending on the config file). It also writes the features included
        in the configuration file, in order to track everything on the same document.
        :return: file
        """
        # Update cases and control lists after preprocessing
        self.update_cohort_lists(self.l_cases, self.output_cases)
        self.update_cohort_lists(self.l_cases_matched, self.output_cases_matched)
        self.update_cohort_lists(self.l_controls, self.output_controls)
        self.update_cohort_lists(self.l_controls_matched, self.output_controls_matched)

        # Write summary of the CAF preprocessing
        json.dump(self.output_log, open(self.output_log_file, 'w'))
        return self.output_log_file

    def update_cohort_lists(self, l_cohort, f_output):
        """
        It updates the final cases and control lists after all pre-processing.
        It will write new self.output_case
        :return:
        """
        with open(f_output, 'w') as f:
            for item in l_cohort:
                if isinstance(item, list):
                    f.write("%s\n" % item[0])
                else:
                    f.write("%s\n" % item)

    def reported_vs_genetic_gender_check(self):
        """
        It filters out samples (in both cohorts) that mismatch the selfreported gender with the
        gender estimated by CGR internal pipeline. This info is in the QC metadata
        (`CGR_predicted_sex` and `SelfDeclaredGender`)
        :return: nothing. updates self.l_controls and self.l_cases (and the log dic)
        """
        # `SelfDeclaredGender` is F or M -- change for Female/Male to compare with predicted gender
        report_gender_cases = ["Female" if x is "F" else "Male" for x in
                               self.metadata_cases['SelfDeclaredGender'].tolist()]
        estimated_gender_cases = self.metadata_cases['CGR_predicted_sex'].tolist()

        mismatch_cases = [i for i, (a, b) in enumerate(zip(report_gender_cases,
                                                           estimated_gender_cases)) if a != b]
        self.output_log['repogen_cases'] = list(np.array(list(self.metadata_cases[self.id_cases]))
                                                [mismatch_cases])
        # Update the list of cases excluding mismatches gender cases
        self.l_cases = np.delete(list(self.metadata_cases[self.id_cases]), mismatch_cases)
        self.metadata_cases = self.metadata_cases.drop(
            self.metadata_cases.index[mismatch_cases]
        )

        report_gender_controls = ["Female" if x is "F" else "Male" for x in
                                  self.metadata_controls['SelfDeclaredGender'].tolist()]
        estimated_gender_controls = self.metadata_controls['CGR_predicted_sex'].tolist()
        mismatch_controls = [(i, j) for i, j in zip(report_gender_controls,
                                                    estimated_gender_controls)
                             if i != j]
        self.output_log['repogen_controls'] = list(np.array(list(self.metadata_controls[
                                                                     self.id_cases]))
                                                   [mismatch_controls])
        # Update the list of controls excluding mismatches gender sequences
        self.l_controls = np.delete(self.metadata_controls[self.id_cases], mismatch_controls)
        self.metadata_controls = self.metadata_controls.drop(
            self.metadata_controls.index[mismatch_controls]
        )

    def checking_matching_features(self, logger):
        """
        Function that calls directly and individually to each matching-feature function
        :param logger: the logger from the main pipeline script
        :return: nothing
        """
        if self.match_gender:
            logger.info("Matching gender within cases and controls...")
            self.matching_gender()
        if self.match_age:
            logger.info("Matching age within cases and controls...")
            self.matching_age()
        if self.match_ancestry:
            logger.info("Matching ancestry within cases and controls...")
            self.matching_ancestry()

    def matching_gender(self):
        """
        Function that finds and filters controls matching the cases and controls as possible
        Using matchControls function from e1071 R library
        https://www.rdocumentation.org/packages/e1071/versions/1.7-1/topics/matchControls
        :return: updated list of cases and controls after being matched between themselves
        """
        # Gender info: used the estimated gender by internal pipeline (in the metadata)
        df_cases = self.metadata_cases[['CGRSequenceID', 'CGR_predicted_sex']]
        df_controls = self.metadata_controls[['CGRSequenceID', 'CGR_predicted_sex']]

        # Selecting controls depending on gender distribution in cases
        matched_controls = epy.phjSelectCaseControlDataset(phjCasesDF=df_cases,
                                                           phjPotentialControlsDF=df_controls,
                                                           phjUniqueIdentifierVarName=
                                                           'CGRSequenceID',
                                                           phjMatchingVariablesList=
                                                           ['CGR_predicted_sex'],
                                                           phjControlsPerCaseInt=
                                                           int(self.cc_ratio),
                                                           phjPrintResults=False)

        # Save matched/selected cases AND controls into matched cases AND controls respectively
        self.l_cases_matched = matched_controls.loc[matched_controls['case'] == 1][[
            'CGRSequenceID']].values.tolist()
        self.l_controls_matched = list(matched_controls.loc[matched_controls['case'] == 0][[
            'CGRSequenceID']].values.tolist())

    def matching_age(self):
        pass

    def matching_ancestry(self):
        pass


