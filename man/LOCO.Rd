% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LOCO.R
\name{LOCO}
\alias{LOCO}
\title{LOCO Class}
\description{
LOCO Class

LOCO Class
}
\details{
Calculates Leave-One-Covariate-Out (LOCO) scores.
}
\examples{
\dontshow{if (requireNamespace("rpart")) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

library(mlr3)

loco = LOCO$new(
  task = tsk("zoo"),
  learner = lrn("classif.rpart"),
  measure = msr("classif.ce")
)

loco$compute()
\dontshow{\}) # examplesIf}
}
\section{Super class}{
\code{\link[xplainvi:FeatureImportanceLearner]{xplainvi::FeatureImportanceLearner}} -> \code{LOCO}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-LOCO-new}{\code{LOCO$new()}}
\item \href{#method-LOCO-compute}{\code{LOCO$compute()}}
\item \href{#method-LOCO-clone}{\code{LOCO$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="xplainvi" data-topic="FeatureImportanceLearner" data-id="print"><a href='../../xplainvi/html/FeatureImportanceLearner.html#method-FeatureImportanceLearner-print'><code>xplainvi::FeatureImportanceLearner$print()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-LOCO-new"></a>}}
\if{latex}{\out{\hypertarget{method-LOCO-new}{}}}
\subsection{Method \code{new()}}{
Creates a new instance of this \link[R6:R6Class]{R6} class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LOCO$new(task, learner, measure, resampling = NULL, features = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{task, learner, measure, resampling, features}}{Passed to \code{FeatureImportanceLearner} for construction.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-LOCO-compute"></a>}}
\if{latex}{\out{\hypertarget{method-LOCO-compute}{}}}
\subsection{Method \code{compute()}}{
A short description...
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LOCO$compute(relation = c("difference", "ratio"))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{relation}}{(character(1)) Calculate \code{"difference"} (default) or \code{"ratio"} of
original scores and scores after permutation}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-LOCO-clone"></a>}}
\if{latex}{\out{\hypertarget{method-LOCO-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LOCO$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
