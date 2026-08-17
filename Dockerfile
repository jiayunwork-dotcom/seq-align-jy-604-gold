FROM golang:1.21-alpine
ENV GOTOOLCHAIN=local CGO_ENABLED=0
WORKDIR /app
COPY go.mod go.sum ./
RUN go mod download
COPY . .
RUN go build -o /app/bin/seq-align .
CMD ["/bin/sh"]
