<template>
  <v-container fluid>
    <v-toolbar dense>
      <v-toolbar-title>
        Edit Group
      </v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn @click="cancel" text color="primary">Cancel</v-btn>
        <v-btn @click="reset" text color="primary">Reset</v-btn>
        <v-btn @click="submit" text :disabled="!valid" color="primary">Save</v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-card class="mt-4 px-4">
      <v-card-text>
        <v-form v-model="valid" ref="form" lazy-validation>
          <v-text-field label="Name" v-model="name" :rules="nameRules" />
          <v-text-field label="Description" v-model="description" />
          <v-text-field label="URL" v-model="url" />
          <v-checkbox label="Open" v-model="isOpen" />
          <v-combobox
            v-model="tags"
            :filter="filter"
            :hide-no-data="!search"
            :items="items"
            :search-input.sync="search"
            hide-selected
            label="Tags"
            multiple
            small-chips
            solo
          >
            <template v-slot:no-data>
              <v-list-item>
                <span class="text-subtitle-1">Create</span>
                <v-chip label small>
                  {{ search }}
                </v-chip>
              </v-list-item>
            </template>
            <template v-slot:selection="{ attrs, item, parent, selected }">
              <v-chip v-if="item === Object(item)" v-bind="attrs" :input-value="selected" label small>
                <span class="pr-2">
                  {{ item.text }}
                </span>
                <v-icon small @click="parent.selectItem(item)">mdi-close</v-icon>
              </v-chip>
            </template>
            <template v-slot:item="{ index, item }">
              <v-text-field
                v-if="editing === item"
                v-model="editing.text"
                autofocus
                flat
                background-color="transparent"
                hide-details
                solo
                @keyup.enter="edit(index, item)"
              ></v-text-field>
              <v-chip v-else dark label small>
                {{ item.text }}
              </v-chip>
              <v-spacer></v-spacer>
              <v-list-item-action @click.stop>
                <v-btn icon @click.stop.prevent="edit(index, item)">
                  <v-icon>{{ editing !== item ? "mdi-pencil" : "mdi-check" }}</v-icon>
                </v-btn>
              </v-list-item-action>
            </template>
          </v-combobox>
        </v-form>
      </v-card-text>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { required } from "@/utils/validators";
import { Component, Vue, Watch } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { IGroupUpdate } from "@/modules/group/models";

@Component
export default class EditGroup extends Vue {
  readonly groupContext = groupModule.context(this.$store);

  readonly nameRules = [required];

  valid = true;
  name = "";
  description: string | null = null;
  url: string | null = null;
  isOpen = false;
  tags: any[] = [];

  editing = null;
  index = -1;
  nonce = 1;
  search = null;

  get items(): any[] {
    const list = this.groupContext.getters.tags;
    return list.map((item) => {
      return {
        text: item,
      };
    });
  }

  get group() {
    return this.groupContext.getters.getGroup(+this.$router.currentRoute.params.groupId);
  }

  async mounted() {
    await Promise.all([
      this.groupContext.actions.getGroup(+this.$router.currentRoute.params.groupId),
      this.groupContext.actions.getTags(),
    ]);
    this.reset();
  }

  reset() {
    this.name = "";
    this.description = "";
    this.url = "";
    this.isOpen = false;
    this.tags = [];
    if (this.$refs.form) {
      (this.$refs.form as any).resetValidation();
    }
    if (this.group) {
      this.name = this.group.name;
      this.description = this.group.description;
      this.url = this.group.url;
      this.isOpen = this.group.is_open;
      this.tags = this.group.tags
        ? this.group.tags.map((item) => {
            return {
              text: item,
            };
          })
        : [];
    }
  }

  cancel() {
    this.$router.back();
  }

  @Watch("tags")
  onTagsChange(val: any[], prev: any[]) {
    if (!val || !prev || val.length === prev.length) {
      return;
    }

    this.tags = val.map((v) => {
      if (typeof v === "string") {
        v = {
          text: v,
        };
        this.items.push(v);
        this.nonce++;
      }
      return v;
    });
  }

  edit(index: number, item) {
    if (!this.editing) {
      this.editing = item;
      this.index = index;
    } else {
      this.editing = null;
      this.index = -1;
    }
  }

  filter(item, queryText: string, itemText: string) {
    const hasValue = (val) => (val != null ? val : "");

    const text = hasValue(itemText);
    const query = hasValue(queryText);

    return text.toString().toLowerCase().indexOf(query.toString().toLowerCase()) > -1;
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const data: IGroupUpdate = {
        name: this.name,
        description: this.description,
        url: this.url,
        is_open: this.isOpen,
        tags: this.tags.map((tag) => tag.text),
      };
      await this.groupContext.actions.updateGroup({
        id: this.group!.id,
        data: data,
      });
      this.$router.back();
    }
  }
}
</script>
