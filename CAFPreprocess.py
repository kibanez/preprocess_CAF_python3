import json
import os
import numpy as np
import pandas as pd
import epydemiology as epy
from datetime import datetime


class CAFPreprocess:

    def __init__(self, hash_cfg):
        # Input parameters
        self.id_cc = hash_cfg.get('id_cc')
        if self.id_cc not in ("CGRSequenceID", "SampleID"):
            raise IOError('The ID used for the list of cases and controls is missing - '
                          'Please specify `CGRSequenceID` or `SampleID`')

        f_cases = hash_cfg.get('list_cases')
        if not os.path.isfile(f_cases):
            raise IOError('The file containing the list of cases %s dose not exist'
                          % f_cases)
        else:
            self.l_cases = self.read_cohort(f_cases, 0)

        f_controls = hash_cfg.get('list_controls')
        if not os.path.isfile(f_controls):
            raise IOError('The file containing the list of controls %s dose not exist'
                          % f_controls)
        else:
            self.l_controls = self.read_cohort(f_controls, 0)

        # matched cases and control lists
        self.l_cases_matched = []
        self.l_controls_matched = []

        f_metadata = hash_cfg.get('metadata_file')
        if not os.path.isfile(f_metadata):
            raise IOError('The QC metadata file %s does not exist' % f_metadata)
        else:
            self.metadata_cases = self.read_metadata(f_metadata,
                                                     self.l_cases,
                                                     self.id_cc)
            self.metadata_controls = self.read_metadata(f_metadata,
                                                        self.l_controls,
                                                        self.id_cc)

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

        self.list_matching = hash_cfg.get('list_matching').split(',')
        if len(self.list_matching) > 0:
            for match_name in self.list_matching:
                if match_name == 'gender':
                    self.match_gender = True
                elif match_name == 'age':
                    self.match_age = True
                elif match_name == "ancestry":
                    self.match_ancestry = True
        else:
            self.match_gender = False
            self.match_age = False
            self.match_ancestry = False

        if self.match_age:
            self.df_cases_age = self.read_cohort_age(f_cases)
            self.df_controls_age = self.read_cohort_age(f_controls)

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

    def read_cohort(self, file_cohort, column):
        """
        Reads a file containing individuals within a cohort (cases/controls) and
        it returns a list of them
        :param file_cohort: path to the file containing the list of individuals
        :param column: the column we want to read
        :return: list of CGRSequenceID
        """
        cohort_data = pd.read_csv(file_cohort, sep=',',
                                  error_bad_lines=False, header=None)
        return cohort_data.iloc[:, int(column)].tolist()

    def read_cohort_age(self, file_cohort):
        """
        Reads a file containing individuals within a cohort (cases/controls) together with the
        age information
        :param file_cohort: path to the file containing the list of individuals + age
        (separated by `,`)
        :return: dataframe containing id and age
        """
        cohort_data = pd.read_csv(file_cohort, sep=',',
                                  error_bad_lines=False, header=None)
        l_ids = cohort_data.iloc[:, 0].tolist()
        l_age = cohort_data.iloc[:, 1].tolist()
        dict_age = {self.id_cc: l_ids, 'age': l_age}
        return pd.DataFrame(dict_age)

    def read_metadata(self, file_metadata, l_cohort, id_selection):
        """
        Reads a file with the QC metadata created by CGR internal bioinformatics pipeline
        :param file_metadata: path to the QC metadata file
        :param l_cohort: list of individuals within that cohort (cases or controls)
        :param id_selection: binary value: "CGRSequenceID" or "SampleID"
        :return: dict with all metadata columns
        """
        l_cohort_str = [str(i) for i in l_cohort]
        qc_metadata = pd.read_csv(file_metadata, sep=',', error_bad_lines=False)
        return qc_metadata.loc[qc_metadata[id_selection].isin(l_cohort_str)]

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
        self.output_log['repogen_cases'] = list(np.array(list(self.metadata_cases[self.id_cc]))
                                                [mismatch_cases])
        # Update the list of cases excluding mismatches gender cases
        self.l_cases = np.delete(list(self.metadata_cases[self.id_cc]), mismatch_cases)
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
                                                                     self.id_cc]))
                                                   [mismatch_controls])
        # Update the list of controls excluding mismatches gender sequences
        self.l_controls = np.delete(self.metadata_controls[self.id_cc], mismatch_controls)
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
        Using epy python library
        https://pypi.org/project/epydemiology/0.1.4/
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
        """
        Function that matches each case in turn and selects the relevant number of control subjects
        from the second dataframe, matching on the list of variables.
        Using epy python library
        https://pypi.org/project/epydemiology/0.1.4/
        :return: updated list of cases and controls after being matched between themselves
        """
        # Selecting controls depending on age distribution in cases
        matched_controls = epy.phjSelectCaseControlDataset(phjCasesDF=
                                                           self.df_cases_age,
                                                           phjPotentialControlsDF=
                                                           self.df_controls_age,
                                                           phjUniqueIdentifierVarName=
                                                           'id',
                                                           phjMatchingVariablesList=
                                                           ['age'],
                                                           phjControlsPerCaseInt=
                                                           int(self.cc_ratio),
                                                           phjPrintResults=False)

        # Save matched/selected cases AND controls into matched cases AND controls respectively
        self.l_cases_matched = matched_controls.loc[matched_controls['case'] == 1][[
            'id']].values.tolist()
        self.l_controls_matched = list(matched_controls.loc[matched_controls['case'] == 0][[
            'id']].values.tolist())

    def matching_ancestry(self):
        """
        Function that matches following c-c ratio cases and controls based on the ancestry
        information
        :return: updated list of cases and controls after being matched between themselves
        """
        # Ancestry info: used the estimated ancestry from peddy from the internal pipeline
        # (in the metadata)
        df_cases = self.metadata_cases[['CGRSequenceID', 'peddy_ancestry_pred']]
        df_controls = self.metadata_controls[['CGRSequenceID', 'peddy_ancestry_pred']]

        # Selecting controls depending on gender distribution in cases
        matched_controls = epy.phjSelectCaseControlDataset(phjCasesDF=df_cases,
                                                           phjPotentialControlsDF=df_controls,
                                                           phjUniqueIdentifierVarName=
                                                           'CGRSequenceID',
                                                           phjMatchingVariablesList=
                                                           ['peddy_ancestry_pred'],
                                                           phjControlsPerCaseInt=
                                                           int(self.cc_ratio),
                                                           phjPrintResults=False)

        # Save matched/selected cases AND controls into matched cases AND controls respectively
        self.l_cases_matched = matched_controls.loc[matched_controls['case'] == 1][[
            'CGRSequenceID']].values.tolist()
        self.l_controls_matched = list(matched_controls.loc[matched_controls['case'] == 0][[
            'CGRSequenceID']].values.tolist())

    def matching_list(self):
        """
        Function that matches a combination of features entered by the user.
        - If the list is empty, only c-c ratio is taken into account and controls are selected
        randomly
        - Features within the list are considered together, controls will be matched depending
        on the values considering all features
        - If the list contains `age` we will do a special analysis on this. since we need to take
        information particularly from the file with the case list
        :return: updated list of cases and controls matched, after being matched between themselves
        """
        # list of variables to match
        if len(self.list_matching) == 0:
            l_matching_variables = None
        else:
            l_matching_variables = self.list_matching
            # Leave `age` at the end of the list, if exists
            if 'age' in l_matching_variables:
                l_matching_variables.remove('age')
                l_matching_variables.append('age')

        matched_controls = epy.phjSelectCaseControlDataset(phjCasesDF=
                                                           self.df_cases_age,
                                                           phjPotentialControlsDF=
                                                           self.df_controls_age,
                                                           phjUniqueIdentifierVarName=
                                                           self.id_cc,
                                                           phjMatchingVariablesList=
                                                           l_matching_variables,
                                                           phjControlsPerCaseInt=
                                                           int(self.cc_ratio),
                                                           phjPrintResults=False)

        # Save matched/selected cases AND controls into matched cases AND controls respectively
        self.l_cases_matched = matched_controls.loc[matched_controls['case'] == 1][[
            self.id_cc]].values.tolist()
        self.l_controls_matched = list(matched_controls.loc[matched_controls['case'] == 0][[
            self.id_cc]].values.tolist())


