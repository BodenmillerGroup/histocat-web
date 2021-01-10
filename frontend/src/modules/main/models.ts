export type ViewMode = "image" | "segmentation" | "data";

export interface AppNotification {
  content: string;
  color?: string;
  showProgress?: boolean;
}
