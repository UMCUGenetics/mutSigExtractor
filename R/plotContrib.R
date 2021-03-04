MUT_SUBTYPE_PARSERS <- list(
   snv=function(x){ gsub('(^\\w\\[)|(\\]\\w$)','',x) },
   indel=function(x){ gsub('[.]\\d+\\+*$','',x) },
   dbs=function(x){ gsub('\\w{2}$','NN',x) },
   sv=function(x){ grep('^[[:upper:]]{3}',x,value=T) }
)

MUT_SUBTYPE_COLORS <- list(
   snv=c(
      "C>A"="#06BAEB",
      "C>G"="#737373",
      "C>T"="#E12825",
      "T>A"="grey",
      "T>C"="#A1CD63",
      "T>G"="#EDC4C5"
   ),

   dbs=c(
      "AC>NN"="#03BCEC",
      "AT>NN"="#0165CB",
      "CC>NN"="#9FCC62",
      "CG>NN"="#006401",
      "CT>NN"="#FE9796",
      "GC>NN"="#E22823",
      "TA>NN"="#FCB066",
      "TC>NN"="#FC8100",
      "TG>NN"="#CC98FD",
      "TT>NN"="#4C0199"
   ),

   indel=c(
      "del.1.C"="#F3BF7B",
      "del.1.T"="#EE8633",
      "ins.1.C"="#B9DA94",
      "ins.1.T"="#569D40",
      "del.2.rep"="#F5CBB7",
      "del.3.rep"="#EE8E72",
      "del.4.rep"="#DE523E",
      "del.5+.rep"="#AD2D25",
      "ins.2.rep"="#D2E0EF",
      "ins.3.rep"="#9DC3DC",
      "ins.4.rep"="#5D97C4",
      "ins.5+.rep"="#2E65A5",
      "del.2.mh"="#E1E1ED",
      "del.3.mh"="#B5B6D5",
      "del.4.mh"="#8585B8",
      "del.5+.mh"="#5D4494",
      "del.mh"="#5D4494"
   ),

   sv=c(
      "DEL"="#9DD1C7",
      "DUP"="#FFFCBB",
      "INV"="#BDBAD7",
      "TRA"="#EB8677"
   )
)

####################################################################################################
#' Plot context or signature contributions
#'
#' @param x A matrix or vector of context or signature contributions
#' @param mode Can be 'signatures' or 'contexts'
#' @param mut.type A string specifying the mutation type Can be 'snv', 'dbs', 'indel'
#' @param group A vector indicating the which group each sample belongs to. If provided, a plot will
#' be produced for each group. If unspecified, the function will assume that all rows of `x` belong
#' to the same group
#'
#' @return A ggplot2 object
#' @export
#'
plotContrib <- function(
   x, mode='contexts', mut.type='auto',
   group=NULL, y.axis.var.scale=F
){

   if(F){
      x=t(sigs.de_novo)
      mode='contexts'
      mut.type='snv'
      group=NULL
      y.axis.var.scale=F
   }

   ## Checks --------------------------------
   require(ggplot2)

   if(length(mut.type)!=1){ stop("`mut.type` must be a string") }
   if(mut.type=='sbs'){ mut.type <- 'snv' }
   if(!(mut.type %in% c('snv','indel','dbs','auto','all'))){
      stop("`mut.type` must be one of the following: 'snv','indel', 'dbs', 'auto','all")
   }

   if(is.vector(x)){ x <- t(x) }
   if(length(colnames(x))==0){ stop("`x` must have colnames if a matrix, or names if a vector") }

   MUT_TYPES_VALID <- list(
      snv=SUBS_CONTEXTS_96,
      dbs=INDEL_CONTEXTS,
      indel=DBS_TYPES$context
   )

   if(mode=='contexts'){

      if(mut.type=='auto'){
         mut_type_tmp <- NULL
         for(i in names(MUT_TYPES_VALID)){
            #i='snv'
            #x=cbind(x, test=NA)
            if(all(MUT_TYPES_VALID[[i]] %in% colnames(x))){
               mut_type_tmp <- i
               break
            }
         }
         if(is.null(mut_type_tmp)){
            stop("`mut.type` is set to 'auto' but could not be inferred from colnames of `x`")
         }
         mut.type <- mut_type_tmp; rm(mut_type_tmp)

      } else {
         if(mut.type=='all'){
            mut_types <- unlist(MUT_TYPES_VALID, use.names=F)
         } else {
            mut_types <- MUT_TYPES_VALID[[mut.type]]
         }

         invalid_mut_types <- colnames(x)[!all(mut_types %in% colnames(x))]
         if(length(invalid_mut_types)>0){
            stop('names of `x` contain invalid or missing mutation types:\n', paste(invalid_mut_types, collapse=', '))
         }
      }
   }

   ## Prep data --------------------------------
   ## Subset for relevant contexts
   df <- as.data.frame(x)
   if(mut.type=='all'){
      df <- df[,unlist(MUT_TYPES_VALID, use.names=F)]
   } else {
      df <- df[, MUT_TYPES_VALID[[mut.type]] ]
   }

   ## Add metadata
   df$sample <- if(length(rownames(x))==0){
      paste0('sample',1:nrow(x))
   } else {
      rownames(x)
   }
   rownames(df) <- NULL

   df$group <- if(!is.null(group)){
      factor(group, unique(group))
   } else {
      'group1'
   }

   df <- reshape2::melt(
      df,
      id.vars=colnames(df)[!(colnames(df) %in% colnames(x))]
   )

   ## Summary stats
   if(nrow(x)!=1){
      df_agg <- aggregate(
         df$value,
         list(group=df$group, context=df$variable),
         function(x){
            c(median=median(x), min=min(x), max=max(x), q1=unname(quantile(x,0.25)), q3=unname(quantile(x,0.75)))
         }
      )
      df_agg <- cbind(
         subset(df_agg, select=-x),
         as.data.frame(df_agg$x)
      )
   } else {
      df_agg <- structure(df, names=c('sample','group','context','median'))
   }
   rm(df)

   ## Get mut subtype (only for contexts) --------------------------------
   if(mode=='contexts'){

      mut_type_lookup <- structure(
         rep(names(MUT_TYPES_VALID), sapply(MUT_TYPES_VALID, length)),
         names=unlist(MUT_TYPES_VALID, use.names=F)
      )

      mut_subtype_lookup <- sapply(names(mut_type_lookup), function(i){
         #i=names(mut_type_lookup)[[1]]
         MUT_SUBTYPE_PARSERS[[ mut_type_lookup[[i]] ]] (i)
      })

      df_agg$context <- as.character(df_agg$context)
      df_agg$mut_type <- unname(mut_type_lookup[df_agg$context])
      df_agg$mut_subtype <- unname(mut_subtype_lookup[df_agg$context])
      df_agg$context <- factor(df_agg$context, unique(df_agg$context))

      ## Merge delmh contexts into one group
      df_agg$mut_subtype_2 <- factor(df_agg$mut_subtype ,unique(df_agg$mut_subtype))
      levels(df_agg$mut_subtype_2)[grep('del.*mh',levels(df_agg$mut_subtype_2))] <- 'del.mh'

      colors <- unlist(unname(MUT_SUBTYPE_COLORS))
      colors <- colors[names(colors) %in% levels(df_agg$mut_subtype_2)]

   } else {
      df_agg$mut_subtype_2 <- factor(df_agg$context, unique(df_agg$context))
      colors <- structure(rep('grey', length(levels(df_agg$mut_subtype_2))))
   }

   ## Plot --------------------------------
   p <- ggplot(df_agg, aes(x=context, y=median))

   ## Facetting
   scales <- 'free_x'
   if(y.axis.var.scale){ scales <- 'free' }

   if(mode=='contexts'){
      if(!is.null(group) & nrow(x)!=1){
         p <- p + facet_grid(group~mut_subtype_2, scales=scales, space='free_x')
      } else {
         p <- p + facet_grid(~mut_subtype_2, scales=scales, space='free_x')
      }
   } else {
      if(!is.null(group) & nrow(x)!=1){
         p <- p + facet_grid(group~., scales=scales, space='free_x')
      }
   }

   p <- p +
      geom_bar(aes(fill=mut_subtype_2),stat='identity') +
      scale_fill_manual(values=colors, name='Mutation type', guide=F)

   ## Error bars
   if(nrow(x)!=1){
      p <- p + geom_linerange(aes(ymin=q1, ymax=q3))
   }

   ## Formatting
   p <- p +
      ylab('Contribution') +
      xlab('Mutation context') +
      theme_bw() +
      theme(
         panel.grid.minor.y=element_blank(),
         #panel.grid.minor.x=element_blank(),
         #panel.grid.major.x=element_blank(),
         axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
         panel.spacing.x=unit(0, "lines"),
         legend.position='left'
      )
   return(p)
}
