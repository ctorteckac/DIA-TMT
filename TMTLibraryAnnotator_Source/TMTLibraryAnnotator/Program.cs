/*
 * Original author: Oliver M. Bernhardt
 * Email: oliver.bernhardt@biognosys.com,
 *
 * Copyright (c) 2020 Oliver Bernhardt
 */

using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TMTLibraryAnnotator {
    //This class takes a Spectronaut spectral library in .xls format and injects TMT reporter ions into the library
    public class Program {

        static void Main(string[] args) {
            string file = args[0];

            System.IO.StreamReader reader = new System.IO.StreamReader(file);
            System.IO.StreamWriter writer = new System.IO.StreamWriter(file.Replace(".xls", "") + "-TMT.tsv");
            string line = reader.ReadLine();
            string[] header = line.Split('\t');
            Dictionary<string, int> headerMap = new Dictionary<string, int>();
            for (int i = 0; i < header.Length; ++i) {
                headerMap[header[i]] = i;
            }

            writer.WriteLine(line);
            System.IO.StreamReader labelsReader = new System.IO.StreamReader("channels.tsv");
            List<double> channels = new List<double>();
            while ((line = labelsReader.ReadLine()) != null) {
                double c;
                if (TrySafeParse(line, false, out c)) {
                    channels.Add(c);
                }
            }

            channels.Sort();

            int sequenceColumn = headerMap["ModifiedPeptide"];
            int midColumn = headerMap["LabeledPeptide"];
            int precZColumn = headerMap["PrecursorCharge"];
            int frgZColumn = headerMap["FragmentCharge"];
            int frgTypeColumn = headerMap["FragmentType"];
            int frgPosColumn = headerMap["FragmentNumber"];
            int q3Column = headerMap["FragmentMz"];
            int excludeColumn = headerMap["ExcludeFromAssay"];

            string lastLineID = "";
            while ((line = reader.ReadLine()) != null) {
                string[] split = line.Split('\t');
                split[sequenceColumn] = split[midColumn];
                string lineID = split[sequenceColumn] + split[precZColumn];

                if (!split[excludeColumn].Equals("TRUE")) {
                    writer.WriteLine(toLine(split));
                }

                if (!lineID.Equals(lastLineID)) {
                    split[frgTypeColumn] = "r";
                    for (int i = 0; i < channels.Count; ++i) {
                        split[frgPosColumn] = (i + 1).ToString();
                        split[frgZColumn] = "1";
                        split[q3Column] = channels[i].ToString();
                        split[excludeColumn] = "FALSE";
                        writer.WriteLine(toLine(split));
                    }
                }

                lastLineID = lineID;
            }

            reader.Close();
            writer.Close();
        }

        private static string toLine(string[] cells) {
            StringBuilder builder = new StringBuilder();
            string sep = "";
            foreach (string c in cells) {
                builder.Append(sep);
                builder.Append(c);
                sep = "\t";
            }
            return builder.ToString();
        }

        #region safe double parsing
        private static char? decimalPointSign = null;
        private static HashSet<string> NAs = new HashSet<string>(new string[] { "NA", "NaN", "N/A", "n. def.", double.NaN.ToString() });
        public static char getSystemDecimalPointSign() {
            return CultureInfo.CurrentCulture.NumberFormat.NumberDecimalSeparator[0];
        }

        private static char getDecimalPointSign(string doubleS) {
            if (doubleS == null) return '.';
            foreach (char c in doubleS) {
                switch (c) {
                    case '.':
                        return '.';
                    case ',':
                        return ',';
                    default:
                        break;
                }
            }
            return '.';
        }

        /// <summary>
        /// This method is meant to be a Culture insensitive way of parsing a string to a double
        /// </summary>
        /// <param name="v"></param>
        /// <param name="setNaAsTrue"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        public static bool TrySafeParse(string v, bool setNaAsTrue, out double result) {
            if (NAs.Contains(v)) {
                result = double.NaN;
                return setNaAsTrue;
            }

            char ds = getDecimalPointSign(v);
            char sds = getSystemDecimalPointSign();
            if (ds != sds) {
                v = v.Replace(ds, sds);
            }
            return double.TryParse(v, out result);
        }
        #endregion
    }
}
