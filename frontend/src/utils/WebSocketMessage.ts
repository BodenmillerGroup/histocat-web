export class WebSocketMessage {
  readonly experimentId: number;
  readonly type: string;
  readonly payload: any;

  constructor(json: object) {
    this.experimentId = json["experiment_id"];
    this.type = json["type"];
    this.payload = json["payload"];
  }
}
