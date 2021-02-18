export class WebSocketMessage {
  readonly projectId: number;
  readonly type: string;
  readonly payload: any;

  constructor(json: any) {
    this.projectId = json["project_id"];
    this.type = json["type"];
    this.payload = json["payload"];
  }
}
