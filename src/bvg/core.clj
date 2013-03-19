(ns bvg.core
 (:use [clojure.java.io :only (copy input-stream file reader writer)])
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

(defn set-up-files
 "Make sure all the necessary files are local for scale."
 [scale]
 (let [args (scale scale-args)
       genes (:genes args)
       expected (:expected args)]
   (get-gzip-file (str genes ".gz") genes)
   (when (= scale :small)
     (get-puzzle-file expected))
   (when (= scale :large)
     (get-gzip-file (str expected ".gz") expected))))

(defn read-lines
 "Return a lazy sequence of lines from file named filename."
 [filename]
 (line-seq (reader filename)))

(defn write-lines
 "Write lines to file."
 [filename lines]
 (with-open [w (writer filename)]
   (doseq [line lines] (.write w (str line "\n")))))

(defn normalized-bit-distance
 "Count of bits by which x and y differ normalized to the expected value."
 [expected x y]
 (let [dist (.bitCount (.xor x y))]
   (if (> dist expected) (- dist expected) (- expected dist))))

(defn make-relate
 "Return a function of two indexes and the bitvectors at those
 indexes that returns the normalized bit distance d between the two
 bitvectors paired with the indexes [d [x y]].  The expected
 mutations per bitvector per generation is the product of the
 bitvector length (scale) and the bitwise probability of mutation."
 [mutation-percentage scale]
 (let [expected (/ (* scale mutation-percentage) 100)
       normalized (partial normalized-bit-distance expected)]
   (fn [x bvx y bvy] (print ".") [(int (normalized bvx bvy)) [x y]])))

(defn find-relations
 "Return a sequence of relations [d [x y]] in population of size
 scale, where d is the normalized-bit-distance between bitvectors x
 and y, and x and y are 0-based line indexes into the genes-file."
 [mutation-percentage scale genes-file]
 (println "find-relations" mutation-percentage scale genes-file)
 (let [relate (make-relate mutation-percentage scale)]
   (letfn [(f-r-s [index s]
             (when-let [f (first s)]
               (print "#")
               (concat
                (map (partial relate index f)
                     (iterate inc (inc index)) (rest s))
                (lazy-seq (f-r-s (inc index) (rest s))))))]
     (f-r-s 0 (map #(BigInteger. % 2) (read-lines genes-file))))))

(defn make-graph
 "Return a graph containing edge [x y] and optionally another graph."
 ([[x y]]
    "A graph with one edge [x y]."
    {x #{y} y #{x}})
 ([graph edge]
    "A new graph composed of graph with edge added to it."
    (println "Adding edge" (inc (count graph)) ":" edge)
    (merge-with into graph (make-graph edge))))

(defn edge-out?
 "True if graph shares exactly one vertex with edge [x y]."
 ([graph [x y]]
    (let [cx (contains? graph x) cy (contains? graph y)]
      (or (and cx (not cy)) (and (not cx) cy)))))

(defn span-population
 "Extract from relations [d [x y]] a connected graph of edges [x y]
 spanning scale vertexes and minimizing the sum of d."
 [scale relations]
 (println "span-population" scale)
 (let [spans? (fn [graph] (>= (count graph) scale))
       comp (fn [[l _] [r _]] (compare l r))
       edges (for [[d e] (sort comp relations)] e)]
   (loop [result (make-graph (first edges)) back (rest edges)]
     (if (spans? result) result
         (let [f (first back) r (rest back)]
           (if (edge-out? result f)
             (recur (make-graph result f) (rest edges))
             (recur result r)))))))

(defn find-leaves
 "Return a sequence of the loneliest vertexes (neighborhood count 1)
  from graph, each paired with its (only) neighbor."
 [graph]
 (for [[v neighbors] graph :when (= 1 (count neighbors))]
   [v (first neighbors)]))

(defn prune-leaves
 "Return a new graph with the vertexes in leaves pruned from graph."
 [leaves graph]
 (println "prune-leaves" (count leaves) (count graph))
 (loop [s (seq leaves) g graph]
   (if-let [[x y] (first s)]
     (recur (rest s) (assoc (dissoc g x) y (disj (g y) x)))
     g)))

(defn extract-genealogy-from-graph
 "Extract a map of child index to parent index from graph.  Mark the
 progenitor's parent as -1."
 [graph]
 (println "extract-genealogy-from-graph" (count graph))
 (loop [result (sorted-map) g graph]
   (let [leaves (find-leaves g)]
     (if (empty? leaves)
       (assoc result (first (first g)) -1)
       (let [new (map (fn [[k v]] (assoc result k v)) leaves)]
         (recur (into result new) (prune-leaves leaves g)))))))

;; This is a hack used to prune relations in the large scale problem.
;; See its conditional use in extract-genealogy-from-graph's related?
;; predicate used to filter the relations sequence.
;;
(defn binomial-standard-deviation
 "The integer closest to standard deviation of binomial (n p),
 which is more than close enough for this problem.
 (standard-deviation 10000 mutation-probability) -=> 40
 (standard-deviation 500 mutation-probability) -=> 8.9 sumpin"
 [n p]
 (Math/round (Math/sqrt (* n p (- 1 p)))))

(defn extract-genealogy-from-genes-file
 "Extract a map of child to parent index from the scale lines of
 bitvector genetic data in genes-file, where mutation-percentage is
 the bitwise probability of mutation and scale is the number of bit
 genes in each bitvector."
 [mutation-percentage scale genes-file]
 (println "extract-genealogy-from-genes-file" genes-file)
 (let [bsd (binomial-standard-deviation scale (/ mutation-percentage 100))
       limit (* 5 bsd)
       related? (if (= scale 500) (constantly true) (fn [[d e]] (< d limit)))]
   (extract-genealogy-from-graph
    (span-population
     scale
     (filter related?
             (find-relations mutation-percentage scale genes-file))))))

(defn write-bitvector-genealogy
 "Read probability from command line along with scale lines of
 genetic data from genes-file, extract a probable genealogy and write
 the result to parents-file."
 [mutation-percentage scale genes-file parents-file]
 (println "write-bitvector-genealogy" parents-file)
 (write-lines
  parents-file
  (for [[child parent]
        (extract-genealogy-from-genes-file
         mutation-percentage scale genes-file)] parent)))

(defn error
 "Throw an Error with message."
 [message]
 (throw (Error. message)))

(defn show-usage
 "Show a usage message for this program."
 []
 (doseq
     [line
      ["Extract the genealogy from a population of ITA bitvectors."
       "Usage: bitvector-genealogy <scale>"
       "Where: <scale> is 'small' or 'large'."
       "Both scales use a 20% bitwise mutation probability."
       "The 'small' scale uses a population of 500 bitvectors,"
       "     each with 500 bits."
       "The 'large' scale uses a population of 10000 bitvectors,"
       "     each with 10000 bits."]]
   (println line)))

(defn -main
 [& args]
 (try (println "-main")
      (when (not= 1 (count args)) (error "Need exactly one argument."))
      (let [scale-kw (keyword (first args))]
        (when (not (and (keyword? scale-kw) (scale-kw scale-args)))
          (error "<scale> must be 'small' or 'large'."))
        (set-up-files scale-kw)
        (let [stuff (scale-kw scale-args)]
          (time (apply write-bitvector-genealogy (:args stuff)))
          (when-let [line (some identity
                                (map (fn [n x y] (when (not= x y) n))
                                     (range)
                                     (read-lines (:expected stuff))
                                     (read-lines (:parents stuff))))]
            (error (str "Results not as expected at line:" line)))))
      (catch Throwable x
        (println x)
        (show-usage))))

;; (-main "small")
;; (-main "large")
