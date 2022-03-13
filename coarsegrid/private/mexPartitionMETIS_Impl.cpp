// Copyright 2020-2022 SINTEF Digital, Mathematics & Cybernetics.
#include <mex.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <metis.h>

namespace {
    using IDX = idx_t;
    using REAL = real_t;
    using IVec = std::vector<IDX>;

    IDX numCells(const mxArray* G)
    {
        const auto* cells = mxGetField(G, 0, "cells");
        return static_cast<IDX>(mxGetScalar(mxGetField(cells, 0, "num")));
    }

    IDX numFaces(const mxArray* G)
    {
        const auto* faces = mxGetField(G, 0, "faces");
        return static_cast<IDX>(mxGetM(mxGetField(faces, 0, "neighbors")));
    }

    IDX numCellFaces(const mxArray* G)
    {
        const auto* cells = mxGetField(G, 0, "cells");
        return static_cast<IDX>(mxGetM(mxGetField(cells, 0, "faces")));
    }

    IVec facepos(const mxArray* G)
    {
        const auto nc = numCells(G);
        auto facepos = IVec(nc + 1);
        const auto* cells = mxGetField(G, 0, "cells");
        const auto* fpos = mxGetField(cells, 0, "facePos");

        if (mxIsInt32(fpos)) {
            const auto* fposVal = static_cast<const int*>(mxGetData(fpos));
            std::transform(fposVal, fposVal + nc + 1, facepos.begin(),
                [](const int pos) -> IDX
            {
                return static_cast<IDX>(pos) - 1;
            });
        }
        else if (mxIsDouble(fpos)) {
            const auto* fposVal = mxGetPr(fpos);
            std::transform(fposVal, fposVal + nc + 1, facepos.begin(),
                [](const double pos) -> IDX
            {
                return static_cast<IDX>(pos) - 1;
            });
        }

        return facepos;
    }

    IVec cellfaces(const mxArray* G)
    {
        const auto ncf = numCellFaces(G);
        auto cf = IVec(ncf);

        const auto* cells = mxGetField(G, 0, "cells");
        const auto* cfaces = mxGetField(cells, 0, "faces");

        if (mxIsInt32(cfaces)) {
            const auto* cfVal = static_cast<const int*>(mxGetData(cfaces));
            std::transform(cfVal, cfVal + ncf, cf.begin(),
                [](const int f)
            {
                return static_cast<IDX>(f) - 1;
            });
        }
        else if (mxIsDouble(cfaces)) {
            const auto* cfVal = mxGetPr(cfaces);
            std::transform(cfVal, cfVal + ncf, cf.begin(),
                [](const double f)
            {
                return static_cast<IDX>(f) - 1;
            });
        }

        return cf;
    }

    IVec neighbours(const mxArray* G)
    {
        const auto nf = numFaces(G);
        const auto ncol = 2;
        auto neigh = IVec(ncol * numFaces(G), -1);

        const auto* faces = mxGetField(G, 0, "faces");
        const auto* N = mxGetField(faces, 0, "neighbors");

        if (mxIsInt32(N)) {
            const auto* NVal = static_cast<const int*>(mxGetData(N));
            for (auto col = 0; col < ncol; ++col) {
                for (auto f = 0*nf; f < nf; ++f) {
                    neigh[2*f + col] = static_cast<IDX>(*NVal++) - 1;
                }
            }
        }
        else if (mxIsDouble(N)) {
            const auto* NVal = mxGetPr(N);
            for (auto col = 0; col < ncol; ++col) {
                for (auto f = 0*nf; f < nf; ++f) {
                    neigh[2*f + col] = static_cast<IDX>(*NVal++) - 1;
                }
            }
        }

        return neigh;
    }

    IVec internalfaces(const IVec& neigh)
    {
        const auto nf = static_cast<IDX>(neigh.size() / 2);
        auto iface = IVec(nf, static_cast<IDX>(-1));

        auto ifaceid = IDX{0};
        for (auto f = 0*nf; f < nf; ++f) {
            if ((neigh[2*f + 0] < 0) ||
                (neigh[2*f + 1] < 0))
            {
                // Boundary face.  Skip.
                continue;
            }

            iface[f] = ifaceid++;
        }

        return iface;
    }

    class Graph
    {
    public:
        explicit Graph(const mxArray* G);

        bool assignWeights(const mxArray* w);

        IDX  numVert() const { return this->numCells_; }
        IDX* vertexPos() { return this->vertexPos_.data(); }
        IDX* vertices() { return this->vertices_.data(); }
        IDX* weights()
        {
            return this->weights_.empty()
                ? nullptr
                : this->weights_.data();
        }

    private:
        struct Edge {
            IDX intface = -1;
            IDX allface = -1;
            IDX cellface = -1;
        };

        IDX numCells_{};
        IDX numIntFaces_{};
        IDX numAllFaces_{};
        IDX numGlobalCellFaces_{};

        IVec vertexPos_{};
        IVec vertices_{};
        IVec weights_{};
        std::vector<Edge> edges_{};

        template <typename IdxFunc>
        void assignWeights(const mxArray* w, IdxFunc&& ifunc);
    };

    Graph::Graph(const mxArray* G)
        : numCells_          (numCells(G))
        , numIntFaces_       (0)
        , numAllFaces_       (numFaces(G))
        , numGlobalCellFaces_(numCellFaces(G))
        , vertexPos_         (numCells_ + 1)
    {
        const auto facePos = facepos(G);
        const auto cf      = cellfaces(G);
        const auto neigh   = neighbours(G);
        const auto iface   = internalfaces(neigh);

        this->vertices_.reserve(facePos.back());
        this->vertexPos_[0] = IDX{0};

        for (auto c = 0*this->numCells_; c < this->numCells_; ++c) {
            for (auto b = facePos[c], e = facePos[c + 1]; b != e; ++b) {
                const auto f       = cf[b];
                const auto ifaceID = iface[f];

                if (ifaceID < IDX{0}) {
                    // Boundary face.  Skip.
                    continue;
                }

                const auto c1    = neigh[2*f + 0];
                const auto other = (c1 != c) ? c1 : neigh[2*f + 1];

                auto edge = Edge{};
                edge.intface  = ifaceID;
                edge.allface  = f;
                edge.cellface = b;

                this->vertices_.push_back(other);
                this->edges_.push_back(edge);

                this->numIntFaces_ =
                    std::max(this->numIntFaces_, ifaceID + 1);
            }

            this->vertexPos_[c + 1] =
                static_cast<IDX>(this->vertices_.size());
        }
    }

    bool Graph::assignWeights(const mxArray* w)
    {
        const auto nweights = static_cast<IDX>(mxGetNumberOfElements(w));
        if (nweights == this->numIntFaces_) {
            this->assignWeights(w, [](const Edge& e) { return e.intface; });
        }
        else if (nweights == this->numAllFaces_) {
            this->assignWeights(w, [](const Edge& e) { return e.allface; });
        }
        else if (nweights == this->numGlobalCellFaces_) {
            this->assignWeights(w, [](const Edge& e) { return e.cellface; });
        }
        else {
            return false;
        }

        return true;
    }

    template <typename IdxFunc>
    void Graph::assignWeights(const mxArray* w, IdxFunc&& ifunc)
    {
        const auto* wPtr = mxGetPr(w);

        auto weights = std::vector<REAL>(this->edges_.size());
        std::transform(this->edges_.begin(), this->edges_.end(),
                       weights.begin(),
            [wPtr, &ifunc](const Edge& e) -> REAL
        {
            return static_cast<REAL>(wPtr[ifunc(e)]);
        });

        auto argmax = std::max_element(weights.begin(), weights.end(),
            [](const double w1, const double w2)
        {
            return std::abs(w1) < std::abs(w2);
        });

        if (argmax == weights.end()) {
            // Empty weights_ array.  Don't do anything else.
            return;
        }

        this->weights_.resize(this->edges_.size());
        const auto maxnorm = std::abs(*argmax);
        std::transform(weights.begin(), weights.end(),
                       this->weights_.begin(),
            [maxnorm](const double x)
        {
            return static_cast<IDX>(10.0e3 * (std::abs(x) / maxnorm));
        });
    }

    class Options
    {
    public:
        Options()
        {
            METIS_SetDefaultOptions(this->array());
            this->objtype() = METIS_OBJTYPE_CUT;
        }

        explicit Options(const mxArray* opt);

        void ufactor(const double value)
        {
            this->options_[METIS_OPTION_UFACTOR] =
                    this->truncate(std::floor(1000.0 * (value - 1.0)));
        }

        void ncuts(const double value)
        {
            this->options_[METIS_OPTION_NCUTS] = this->truncate(value);
        }

        void seed(const double value)
        {
            this->options_[METIS_OPTION_SEED] = this->truncate(value);
        }

        void no2hop(const double value)
        {
            this->options_[METIS_OPTION_NO2HOP] = this->truncate(value);
        }

        void minconn(const double value)
        {
            this->options_[METIS_OPTION_MINCONN] = this->truncate(value);
        }

        IDX* array()
        {
            return this->options_.data();
        }

    private:
        using option = void (Options::*)(const double value);

        std::array<IDX, METIS_NOPTIONS> options_{};
        static std::map<std::string, option> select_;

        IDX& objtype()
        { return this->options_[METIS_OPTION_OBJTYPE]; }

        void assign(const char* optname, const mxArray* value);
        double optionValue(const mxArray* value);
        IDX truncate(const double value) const;
    };

    Options::Options(const mxArray* opt)
       : Options{}
    {
        const auto nfld = mxGetNumberOfFields(opt);
        for (auto fld = 0*nfld; fld < nfld; ++fld) {
            const auto* optname = mxGetFieldNameByNumber(opt, fld);
            this->assign(optname, mxGetFieldByNumber(opt, 0, fld));
        }
    }

    void Options::assign(const char* optname, const mxArray* value)
    {
        auto opt = Options::select_.find(optname);
        if (opt == Options::select_.end()) {
            // No such option.
            return;
        }

        (this->*opt->second)(this->optionValue(value));
    }

    double Options::optionValue(const mxArray* value)
    {
        return mxGetScalar(value);
    }

    IDX Options::truncate(const double value) const
    {
        return static_cast<IDX>(value);
    }

    std::map<std::string, Options::option> Options::select_
    {
        { "ufactor", &Options::ufactor },
        { "ncuts",   &Options::ncuts   },
        { "seed",    &Options::seed    },
        { "no2ho",   &Options::no2hop  },
        { "minconn", &Options::minconn },
    };

    IDX num_blocks(const mxArray* nblk)
    {
        return static_cast<IDX>(mxGetScalar(nblk));
    }

    namespace verify_input {
        int verify_faces_structure(mxArray *faces)
        {
            // Shallow structural inspection only.  Assume valid fields...
            return mxGetFieldNumber(faces, "neighbors") >= 0;
        }

        bool verify_cells_structure(mxArray *cells)
        {
           // Shallow structural inspection only.  Assume valid fields...
           bool ok = mxGetFieldNumber(cells, "facePos"  ) >= 0;

           ok = ok && (mxGetFieldNumber(cells, "faces"  ) >= 0);

           return ok;
        }

        bool is_grid(const mxArray *G)
        {
            bool nodes_ok = false, faces_ok = false, cells_ok = false;
            int field_no;

            mxArray *pm;
            if (mxIsStruct(G)) {
                nodes_ok = mxGetFieldNumber(G, "nodes") >= 0;

                if ((field_no = mxGetFieldNumber(G, "faces")) >= 0) {
                    pm = mxGetFieldByNumber(G, 0, field_no);
                    faces_ok = verify_faces_structure(pm);
                }

                if ((field_no = mxGetFieldNumber(G, "cells")) >= 0) {
                    pm = mxGetFieldByNumber(G, 0, field_no);
                    cells_ok = verify_cells_structure(pm);
                }
            }

            return nodes_ok && faces_ok && cells_ok;
        }

        bool is_num_blocks(const mxArray* nblk)
        {
            return mxIsNumeric(nblk) && (mxGetNumberOfElements(nblk) == 1);
        }

        bool is_options(const mxArray* opt)
        {
            return mxIsStruct(opt) && (mxGetNumberOfElements(opt) == 1);
        }

        bool is_weighting_array(const mxArray* w)
        {
            return mxIsNumeric(w) && !mxIsEmpty(w);
        }

        bool args(int nlhs, int nrhs, const mxArray *prhs[])
        {
             auto ok = (nrhs > 1) && (nrhs < 5) && (nlhs < 2)
                       && is_grid(prhs[0]) && is_num_blocks(prhs[1]);

             if (ok && (nrhs > 2)) {
                 ok = is_options(prhs[2]);
             }

             if (ok && (nrhs > 3)) {
                 ok = is_weighting_array(prhs[3]);
             }

             return ok;
        }
    } // verify_input

    void reportErrors(const int nlhs, const int nrhs, const mxArray* prhs[])
    {
        std::ostringstream msg;
        msg << "Invalid MEX Function Call\n";
        if (nrhs < 2) {
            msg << "  - Too few input arguments\n"
                << "    Must provide at least a graph ('G') and "
                << "a number of coarse blocks ('numBlocks')";
        }
        else if (nrhs > 4) {
            msg << "  - Too many input arguments\n"
                << "    Must provide at most a graph ('G'), a "
                << "number of coarse blocks ('numBlocks'), an "
                << "options structure ('opt'), and a weighting "
                << "array ('w')";
        }

        if (nlhs > 1) {
            msg << "\n  - Unsupported number of return values "
                << nlhs << "\n    Must have exactly one return value";
        }

        if (! verify_input::is_grid(prhs[0])) {
            msg << "\n  - Input graph ('G') is not a valid graph";
        }

        if (! verify_input::is_num_blocks(prhs[1])) {
            msg << "\n  - Input number of blocks ('numBlocks') is not a valid scalar";
        }

        if ((nrhs > 2) && ! verify_input::is_options(prhs[2])) {
            msg << "\n  - Input options structure ('opt') is not valid";
        }

        if ((nrhs > 3) && ! verify_input::is_weighting_array(prhs[3])) {
            msg << "\n  - Input weighting array ('w') is not valid";
        }

        mexErrMsgTxt(msg.str().c_str());
    }

    mxArray* makePartitionVector(const IVec& p)
    {
        auto* partition = mxCreateDoubleMatrix(p.size(), 1, mxREAL);

        std::transform(p.begin(), p.end(), mxGetPr(partition),
            [](const IDX blockID) -> double
        {
            return blockID + 1;
        });

        return partition;
    }

    namespace partition_metis
    {
        IVec KWay(Graph& graph, IDX numBlocks, Options& options)
        {
            auto nvtx = graph.numVert();
            auto partition = IVec(nvtx);

            IDX   ncon    = 1;
            IDX*  vwgt    = nullptr;  // No vertex weights
            IDX*  vsize   = nullptr;  // No vertex size
            REAL* tpwgts  = nullptr;
            REAL* ubvec   = nullptr;
            IDX   edgecut = 0;

            const int ret =
                METIS_PartGraphKway(&nvtx, &ncon, graph.vertexPos(),
                                    graph.vertices(), vwgt, vsize,
                                    graph.weights(), &numBlocks, tpwgts,
                                    ubvec, options.array(), &edgecut,
                                    partition.data());

            if (ret == METIS_OK) {
                return partition;
            }
            else if (ret == METIS_ERROR_INPUT) {
                std::string msg = "METIS_PartGraphKway(): Input error";
                throw std::domain_error(msg);
            }
            else if (ret == METIS_ERROR_MEMORY) {
                std::string msg = "METIS_PartGraphKway(): Allocation Failure";
                throw std::runtime_error(msg);
            }
            else {
                std::string msg = "METIS_PartGraphKway(): Other Error";
                throw std::runtime_error(msg);
            }
        }
    }
} // Anonymous

// p = mexPartitionMETIS(G, numBlocks, opt, w)
void mexFunction(int nlhs, mxArray*       plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (verify_input::args(nlhs, nrhs, prhs)) {
       auto graph = Graph{ prhs[0] };
       const auto numBlocks = num_blocks(prhs[1]);
       auto opt = (nrhs > 1) ? Options { prhs[2] } : Options{};
       if (nrhs > 3) {
           graph.assignWeights(prhs[3]);
       }

       const auto p = partition_metis::KWay(graph, numBlocks, opt);
       plhs[0] = makePartitionVector(p);
    }
    else {
        reportErrors(nlhs, nrhs, prhs);
    }
}
