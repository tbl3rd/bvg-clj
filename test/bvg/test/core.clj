(ns bvg.test.core
  (:use [bvg.core])
  (:use [clojure.test])
  (:use [clojure.java.io])
  (:import [java.util.zip GZIPInputStream]))

(defn puzzle-file
  "Return the URL of the puzzle file named leafname."
  [leafname]
  (str "http://www.itasoftware.com/careers/work-at-ita/PuzzleFiles/" leafname))

(defn add-args
  "Add test boilerplate to stuff."
  [stuff]
  (let [more (conj stuff [:percent 20])
        args (into [] (map more [:percent :scale :genes :parents]))]
    (conj more [:args args])))

(def scale-args
  {:small (add-args {:scale 500
                     :genes "bitvectors-genes.data.small"
                     :parents "bitvectors-parents.data.small"
                     :expected "bitvectors-parents.data.small.txt"})
   :large (add-args {:scale 10000
                     :genes "bitvectors-genes.data"
                     :parents "bitvectors-parents.data"
                     :expected "bitvectors-parents.data.txt"})})

(defn get-puzzle-file
  "Get the puzzle file named filename if necessary."
  [filename]
  (let [out-file (file filename)]
    (when (not (.exists out-file))
      (copy (input-stream (puzzle-file filename)) out-file)
      (println "Got:" filename))))

(defn get-gzip-file
  "Get the gzipped file in-name if necessary and uncompress it into
  the file named out-name if necessary."
  [in-name out-name]
  (let [out-file (file out-name) in-file (file in-name)]
    (when (not (.exists out-file))
      (if (not (.exists in-file)) (get-puzzle-file in-name))
      (let [in-stream (GZIPInputStream. (input-stream in-file))]
        (copy in-stream out-file)
        (println "Got:" out-name)))))

(defn test-stuff
  "True if testing against stuff succeeds.  Otherwise false."
  [stuff]
  (println (:args stuff))
  (time (apply write-bitvector-genealogy (:args stuff)))
  (not-any? true?
            (map not=
                 (read-lines (:expected stuff))
                 (read-lines (:parents stuff)))))

(defn setup-puzzle-files
  []
  (let [small-genes (:genes small-test-stuff)
        small-expected (:expected small-test-stuff)
        large-genes (:genes large-test-stuff)
        large-expected (:expected large-test-stuff)]
    (get-gzip-file (str small-genes ".gz") small-genes)
    (get-gzip-file (str large-genes ".gz") large-genes)
    (get-puzzle-file small-expected)))

(defn cleanup-puzzle-files
  []
  (let [delete (fn [fname]
                 (delete-file fname :stoically)
                 (println "Deleted:" fname))
        filenames
        [(:genes small-test-stuff)
         (:parents small-test-stuff)
         (:expected small-test-stuff)
         (:genes large-test-stuff)
         ; (:expected large-test-stuff)   ; computed in C
         (str (:genes small-test-stuff) ".gz")
         (str (:expected small-test-stuff) ".gz")
         (str (:genes large-test-stuff) ".gz")]]
    (doseq [fname filenames] (delete fname))))

(defn small-scale
  []
  (println "Running small test.")
  (is (test-stuff small-test-stuff))
  (println "Small test succeeded."))

(defn large-scale
  []
  (println "Running large test.")
  (is (test-stuff large-test-stuff))
  (println "Large test succeeded."))

(deftest large-then-small
  (setup-puzzle-files)
  (large-scale)
  (small-scale)
  (cleanup-puzzle-files))

;; (setup-puzzle-files)
;; (cleanup-puzzle-files)
