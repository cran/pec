d <- SimSurv(100,cens=FALSE)
d <- SimSurv(1000)
cfit <- coxph(Surv(time,status)~X1+X2,data=d)
a <- pseudoPec(list(cfit),formula=Surv(time,status)~1,data=d)
b <- pec(list(cfit),formula=Surv(time,status)~1,data=d)

plot(a)
plot(b,add=T,lty=2)

data(GBSG2)
model2=coxph(Surv(time,cens)~age+tsize+tgrade+pnodes+progrec+estrec+menostat,data=GBSG2)
a <- pseudoPec(list(model2),formula=Surv(time,cens)~1,data=GBSG2)
b <- pec(list(model2),formula=Surv(time,cens)~1,data=GBSG2)

pdf("~/tmp/test.pdf")
plot(a)
plot(b,add=T,lty=2)
dev.off()

